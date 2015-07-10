#!/usr/bin/python

import sys
import os.path
import argparse
import ConfigParser
import json
import re
import smtplib
from datetime import date, datetime

# Script to report the datum status as determined by the CORS_Analyst (analyse_cors) script
# Uses a configuration file defining the options to use for sending information

class Mailer( object ):

    def __init__( self, config ):
        server=config('smtp_server')
        port=config('smtp_port')
        user=''
        pwd=''
        smtp_file=config('smtp_config_file')
        if smtp_file:
            with open(smtp_file) as smtpf:
                for l in smtpf:
                    parts=l.strip().split(' ',1)
                    if len(parts) == 2:
                        item, value=parts
                        if item == 'server': server=value
                        if item == 'port': port=value
                        if item == 'user': user=value
                        if item == 'password': pwd=value
        self._server=server
        self._port=port
        self._user=user
        self._pwd=pwd

    def send( self, address_from, address_to, subject, message ):
        receivers=address_to.split()
        message=("From: "+address_from+
                 "\nTo: "+address_to+
                 "\nSubject: "+subject+
                 "\n\n"+message)
        try:
            port = int(self._port) if self._port else None
            smtp=smtplib.SMTP(self._server, port)
            if self._user:
                smtp.login(self._user,self._pwd)
            smtp.sendmail(address_from,receivers,message)
        except:
            msg=str(sys.exc_info()[1])
            raise RuntimeError("Unable to send message: "+subject+"\n"+msg)


def get_configuration( config_file, config_section ):
    '''
    Return a function to get configuration items
    '''
    if not config_file:
        config_file=os.path.basename(__file__)
        config_file=os.path.splitext(config_file)[0]+'.cfg'
    cfg=ConfigParser.ConfigParser()
    cfg.read(config_file)
    def getConfigItem( item, default=None, lookup=None ):
        lookup = lookup or {};
        try:
            value=cfg.get(config_section,item)
            if value == '':
                value=default
            elif value != '':
                # Replace lines with just a . with a blank line
                value=re.sub(r'\n\.\n',"\n\n",value,re.S)
                value=re.sub(r'\n\.\n',"\n\n",value,re.S)
                value=re.sub(r'^.\n',"\n",value,re.S)
                value=re.sub(r'\n.$',"\n",value,re.S)
                # Replace {xxx} interpolated strings with lookup, environment, config items
                value=re.sub(r'\{(\w+)\}',lambda m: str(lookup.get(m.group(1),m.group(0))),value)
                value=re.sub(r'\{(\w+)\}',lambda m: os.environ.get(m.group(1),m.group(0)),value)
                value=re.sub(r'\{(\w+)\}',lambda m: cfg.get(config_section,m.group(1),m.group(0)),value)
        except:
            value=default
        return value
    return getConfigItem

def loadDatumStatus( statusFile ):
    with open(statusFile) as sf:
        status=json.load(sf)
        return status

def getWarnings( status ):
    tests=['gdb_offset','scm_offset']
    warnings=[]
    for station in status["station_summary"]:
        for test in tests:
            if station[test]['status_message'] is not None:
                warning=station[test]
                warning['warning_type']=test
                warning['code']=station['code']
                warning['offset_test_date']=station['offset_test_date']
                warnings.append(warning)
    return warnings
                
def filterNewWarnings( warnings,notifiedFile,reminderDays,expiryDays):
    today=date.today()
    dateformat="%Y-%m-%d"
    todaystr=today.strftime(dateformat)
    notified=[]
    try:
        with open(notifiedFile) as nf:
            notified=json.load(nf)
    except:
        pass

    # Filter out already notified warnings, and update notification status
    new_warnings=[]
    for w in warnings:
        try:
            match=(w['code'],w['warning_type'])
            found=filter(lambda x: (x['code'],x['warning_type'])==match,notified)
            found=found[0] if len(found)>0 else None
            if found:
                found['last_seen']=todaystr
                notedate=datetime.strptime(found['last_notified'],dateformat).date()
                noteage=(today-notedate).days
                if noteage > reminderDays:
                    w['status_message']='Reminder: '+w['status_message']
                    new_warnings.append(w)
                    found['last_notified']=todaystr
            else:
                new_warnings.append(w)
                notified.append({
                    'code': w['code'],
                    'warning_type': w['warning_type'],
                    'last_notified': todaystr,
                    'last_seen': todaystr
                    });
        except:
            pass

    # Filter out expired notes
    current_notes=[]
    for n in notified:
        try:
            lastdate=datetime.strptime(n['last_seen'],dateformat).date()
            age=(today-lastdate).days
            if age < expiryDays:
                current_notes.append(n)
        except:
            pass
    notified=current_notes

    return new_warnings, notified

def saveNotificationStatus( notifiedFile, notified ):
    if not notifiedFile:
        return
    try:
        with open(notifiedFile,"w") as nf:
            nf.write(json.dumps(notified ,sort_keys=True,indent=2))
    except:
        pass


def main():
    parser=argparse.ArgumentParser('Report on CORS datum integrity status')
    parser.add_argument('config_file',nargs='?',help="Configuration file for script")
    parser.add_argument('-s','--config-section',default='daily_report',help="Configuration section to use")
    parser.add_argument('-c','--dump-config-file',action='store_true',help="Print an example configuration file")
    args=parser.parse_args()

    if args.dump_config_file:
        cfgfile=os.path.splitext(__file__)[0]+'.cfg'
        with open(cfgfile) as cfgf:
            print cfgf.read()
            sys.exit()

    config=get_configuration(args.config_file,args.config_section)


    datumStatus=loadDatumStatus(config('cors_datum_status_file'))
    warnings=getWarnings(datumStatus)

    # If we are looking for new warnings, then check with 
    # notified warnings file

    notified=None
    notificationFile=config('notified_warnings_file')
    if len(warnings) > 0:
        if notificationFile is not None:
            reminderDays=int(config('notification_reminder_frequency_days',default='30'))
            expiryDays=int(config('notification_expiry_days',default='30'))
            warnings,notified=filterNewWarnings(warnings,notificationFile,reminderDays,expiryDays)

    warning_message=config('no_warning_message')
    if len(warnings) > 0:
        warning_message="\n".join([config('warning_message',lookup=w) for w in warnings])

    if not warning_message:
        sys.exit()

    datumStatus['warnings']=warning_message
    datumStatus['nwarnings']=str(len(warnings))
    message=config('message',lookup=datumStatus)
    title=config('title',lookup=datumStatus)

    from_address=config('from_address')
    to_address=config('to_address')

    try:
        mailer=Mailer(config)
        mailer.send(from_address, to_address, title, message )
        if notified is not None and notificationFile is not None:
            saveNotificationStatus(notificationFile,notified)
    except:
        print str(sys.exc_info()[1])


if __name__=="__main__":
    main()

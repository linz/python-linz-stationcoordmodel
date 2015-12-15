#!/usr/bin/python

import argparse
import re
import os
import os.path
import md5
import sys
import zipfile

def main():
    parser=argparse.ArgumentParser(description="Create a PositioNZ-PP station coordinate model upload file")
    parser.add_argument('zip_file',help="Name of zip file to create")
    parser.add_argument('spm_xml_file',nargs='*',help='Station prediction model XML files to upload CCCC.xml')
    parser.add_argument('-s','--stn-dir',help='Directory containing station coordinate XML files')
    parser.add_argument('-r','--remove-file',nargs='*', help='The station code of files to remove from positionzpp')
    parser.add_argument('-o','--overwrite',action='store_true', help='Allow overwriting an existing zip file')
    parser.add_argument('-v','--verbose',action='store_true', help='Be verbose!')

    args=parser.parse_args()
    verbose=args.verbose

    upload_key='GNSSSPM'
    filere=re.compile(r'^(\w{4})\.xml$')

    if os.path.exists(args.zip_file) and not args.overwrite:
        print("Use -o to overwrite existing zip file "+args.zip_file)
        sys.exit()

    codes=[]
    files=[]
    spmhash=''
    spmdir=args.stn_dir or '.'


    xmlfiles=list(args.spm_xml_file)
    if len(xmlfiles) == 0:
        xmlfiles=[f for f in os.listdir(spmdir) if filere.match(f)]

    for sfile in xmlfiles:
        if not sfile.endswith('.xml'):
            sfile=sfile+'.xml'
        source=os.path.join(spmdir,sfile)
        if not os.path.exists(source):
            source=sfile
        if not os.path.exists(source):
            print(sfile,'does not exist!')
            continue
        filename=os.path.basename(source)
        match=filere.match(filename)
        if not match:
            print(sfile,' is not a valid filename for a station prediction model')
            continue
        code=match.group(1).upper()
        if verbose:
            print('Adding model for {0} from {1}'.format(code,source))
        codes.append(code)
        with open(source) as sf:
            data=sf.read()
            m=md5.new()
            m.update(upload_key)
            m.update(data)
            spmhash=spmhash+code+".xml "+m.hexdigest()+"\n"
        files.append({'file':code+'.xml','data':data})

    remove_files=args.remove_file or []
    for rcode in remove_files:
        if not re.match(r'^\w{4}$',rcode):
            print("Invalid code".rcode," for removal")
            continue
        rcode=rcode.upper()
        if rcode in codes:
            print("Cannot remove {0} as already used".format(rcode))
            continue
        m=md5.new()
        m.update(upload_key)
        m.update("REMOVE")
        m.update(rcode)
        m.update('.xml')
        spmhash=spmhash+rcode+".xml "+m.hexdigest()+"\n"
        if verbose:
            print("Removing file "+rcode)

    zipfilename=args.zip_file
    if not zipfilename.endswith('.zip'):
        zipfilename=zipfilename+'.zip'
    zip=zipfile.ZipFile(zipfilename,'w',zipfile.ZIP_DEFLATED)
    zip.writestr('spm.hash',spmhash)
    for f in files:
        zip.writestr(f['file'],f['data'])
    zip.close()
    if verbose:
        print('Zip file {0} ready for upload to PositioNZ-PP'.format(zipfilename))

if __name__ == '__main__':
    main()

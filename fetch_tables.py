#!/usr/bin/env python3

## Fetch the CERNbox tables and save them in the corresponding directories.

import subprocess
import xml.etree.ElementTree as ET
import sys
import os
import argparse

paths = { # target dir, CERNbox hash
    'mesons' :  ['ALP_production/tab_mesons', 'YBs0coZnaFhp6z5'],
    'prod' :  ['tab_prod', 'MhIAbEQ7GlznA3j'],
    'decay' :  ['tab_decay', '53FYA4vILiEdi9c'],
    'toPlot' :  ['tab_toPlot', 'ic6Siwyl7u0ijBC'],
}

def fetch_folder(hash, replace_files, out_path, in_path = '/'):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    url = 'https://cernbox.cern.ch/remote.php/dav/public-files/' + hash + in_path

    # Call curl to get the XML response
    curl_command = ["curl", "-X", "PROPFIND", url]
    response = subprocess.run(curl_command, capture_output=True, text=True)

    # Parse the XML response and look for subdirectories
    root = ET.fromstring(response.stdout)
    subdir_paths = root.findall('.//d:response/d:href', {'d': 'DAV:'})

    # Extract and print subdirectory names
    for paths in subdir_paths:
        subdir = paths.text.split(hash)[-1]
        if subdir == in_path: continue
        if subdir[-1] == '/':
            if not os.path.exists(out_path+subdir):
                os.makedirs(out_path+subdir)
            fetch_folder(hash, replace_files, out_path, in_path = subdir)
        else:
            if not os.path.exists(out_path + subdir):
                print(f"Downloading {out_path + subdir}")
                download_command = ["curl", "-s", "-o", out_path + subdir, url + subdir.split('/')[-1]]
                subprocess.run(download_command)
            elif replace_files == True:
                print(f"Replacing {out_path + subdir}")
                download_command = ["curl", "-s", "-o", out_path + subdir, url + subdir.split('/')[-1]]
                subprocess.run(download_command)
            else:
                reading = input("File already exists: " + out_path + subdir + ". Replace? (y/n)")
                if reading == 'y':
                    download_command = ["curl", "-s", "-o", out_path + subdir, url + subdir.split('/')[-1]]
                    subprocess.run(download_command)
                elif reading == 'n':
                    continue
                else:
                    print(f"Invalid input. ({reading}) Exitting.")
                    exit()


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)
        
    tabs = ['all']
    tabs.extend(list(paths.keys()))

    parser = argparse.ArgumentParser(description='Fetch the CERNbox tables and save them in the corresponding directories.')
    parser.add_argument("-tab","--tables",  default="", type=str, nargs='*', help="Tables (case sensitive): "+", ".join(tabs))
    parser.add_argument("-f","--replace", action='store_true', help="Replace existing files (default false).")
    parser.set_defaults(replace=False)

    args = parser.parse_args()

    if not args.tables:
        print(f"Please select tables to fetch. (Options: {', '.join(tabs)}) Exitting.")
        exit()

    if 'all' in args.tables:
        if not args.tables == ['all']:
            print(f"Either select \'all\' tables or put a list of individual tables. Exitting.")
            exit()

        for key in paths.keys():
            fetch_folder(paths[key][1], args.replace, paths[key][0], in_path = '/')
    else:
        for tables in args.tables:
            if tables not in tabs:
                print(f"Invalid table name. ({tables}) Exitting.")
                exit()
            fetch_folder(paths[tables][1], args.replace, paths[tables][0], in_path = '/')
            
if __name__ == "__main__":
    sys.exit(main())
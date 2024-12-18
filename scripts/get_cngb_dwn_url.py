#!/share/data1/software/miniconda3/envs/jinxin/bin/python3
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-07-13, 12:18:49
# modified date: 2024-03-30, 12:18:49

'''
Jinxin Meng, 20240330
Generate the run url under the CNGBdb proj.
'''

import re, sys, time
import urllib.request

if len(sys.argv) != 3:
    print("Usage: get_cngb_dwn_url.py [proj_url] [out_f]")
    print(" [proj_url] such as https://ftp.cngb.org/pub/CNSA/data5/CNP0002106/")
    sys.exit()

def main(cnp_url, out_f, sec):
    out_f = open(out_f, "w")
    out_f.write("proj\tcns\tcnx\tcnr\turl\n")

    if cnp_url[-1] != "/":
        cnp_url = cnp_url + "/"
    
    cnp_url_open = urllib.request.urlopen(cnp_url)
    time.sleep(sec)
    cnp_url_html = cnp_url_open.read().decode()
    cns_list = re.findall(r'href="(CNS\d+)/"', cnp_url_html)
    proj_name = re.findall('(CNP\d+)', cnp_url)[0]

    for cns in cns_list:
        cns_url = cnp_url + cns + "/"
        cns_url_open = urllib.request.urlopen(cns_url)
        time.sleep(sec)
        cns_url_html = cns_url_open.read().decode()
        cnx_list = re.findall(r'href="(CNX\d+)/"', cns_url_html)
        
        for cnx in cnx_list:
            cnx_url = cns_url + cnx + "/"
            cnx_url_open = urllib.request.urlopen(cnx_url)
            time.sleep(sec)
            cnx_url_html = cnx_url_open.read().decode()
            cnr_list = re.findall(r'href="(CNR\d+)/"', cnx_url_html)
            
            for cnr in cnr_list:
                cnr_url = cnx_url + cnr + "/"
                cnr_url_open = urllib.request.urlopen(cnr_url)
                time.sleep(sec)
                cnr_url_html = cnr_url_open.read().decode()
                file_list = re.findall(r'href="(.*f[a|q].*)"\srel', cnr_url_html)

                for f in file_list:
                    file_url = cnr_url + f
                    out_f.write("\t".join([proj_name, cns, cnx, cnr, file_url]) + "\n")
    out_f.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sec = 3)

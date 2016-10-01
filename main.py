#!/usr/bin/python
# -*- coding: utf-8 -*-

# import oauth2 as oauth
# import urllib2 as urllib

from pln.pln import Pln

if __name__ == '__main__':
    pln = Pln()
    pln.read_config_json_file()
    pln.read_input_file()
    pln.extract_motifs(pln._motifs_data)
    pln.call_psimod()
    pln.call_prosite()
    pln.print_to_file(pln.motif_and_modification_list, pln.prosite_result, pln.psimod_result)


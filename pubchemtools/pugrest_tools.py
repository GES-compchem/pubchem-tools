#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:09:49 2021

@author: mattia
"""

from dataclasses import dataclass
import requests
from requests.exceptions import HTTPError
import time


# %% debug HTTP performance
# import logging
# logging.basicConfig(level=logging.DEBUG)

# %% Reference dataclass definition
@dataclass
class References:
    reference_id: str = None
    source_name: str = None
    hazard_codes: list = None
    companies: int = None
    notifications: int = None


# %% HTTP Catch errors and reuse session in HTTP request
def HTTP_request(url, session, data=None, other=None):
    if session is None:
        session = requests.Session()

    for i in range(1,4):
        try:
            response = session.post(url, data=data)

            # If the response was successful, no Exception will be raised
            response.raise_for_status()

            if ('-' not in str(other)) and (other is not None):
                print(f"\x1b[1;31;80mData found for {other}\x1b[1;31;0m")

        except HTTPError as http_err:
            if "PUGVIEW.NotFound" in str(http_err):
                print(f"No GHS data found for compound CID {other}")
            if "PUGREST.NotFound" in str(http_err):
                print(f"Compound {other} not found in PUBCHEM")
            break
            #print(f"HTTP error occurred: {http_err}")  # Python 3.6
        except Exception as err:
            #print(f"Other error occurred: {err}")  # Python 3.6
            print(f"\n Connection error. Retrying {other} - Attempt {i}/3")
            time.sleep(1)
        else:
            return response
    else:
        print(f'Three attempts have failed for {other}')


# %%
def cas_to_cid(cas_number, session=None):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/cids/JSON"

    response = HTTP_request(url, session)

    response_json = response.json()

    return response_json["IdentifierList"]["CID"]


def inchi_to_cid(inchi, session=None):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON"

    inchi = "inchi=" + inchi

    response = HTTP_request(url, session, data=inchi)

    response_json = None

    try:
        response_json = response.json()
    except:
        return [None]

    return response_json["IdentifierList"]["CID"]


def inchikey_to_cid(inchikey, session=None):

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"

    start_time = time.time()

    response = HTTP_request(url, session, other=inchikey)

    response_json = None

    try:
        response_json = response.json()
    except:
        return [None]

    return response_json["IdentifierList"]["CID"]


# %%
def cid_to_ghs(cid_number, session=None):

    if cid_number == None:
        #return {None: "Compound not found"}
        return {-1: References(reference_id=-1, source_name="Compound not found")}

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid_number}/JSON/?heading=GHS+Classification"

    response = HTTP_request(url, session, other=cid_number)

    try:
        response_json = response.json()
    except:
        #return {None: "No hazard phrases found"}
        return {-2: References(reference_id=-2, source_name="No hazard phrases found")}

    # data from Pubchem
    ref_sources = response_json["Record"]["Reference"]
    ref_data = response_json["Record"]["Section"][0]["Section"][0]["Section"][0][
        "Information"
    ]

    references = {}

    for source in ref_sources:
        ref_id = source["ReferenceNumber"]
        source_name = source["SourceName"]
        references[ref_id] = References(reference_id=ref_id, source_name=source_name)

    for ref in ref_data:
        ref_id = ref["ReferenceNumber"]
        ref_name = ref["Name"]

        if ref_name == "GHS Hazard Statements":
            #print('\n \x1b[1;34;80mGHS Hazard Statements\x1b[1;34;0m \n')
            hazard_codes = []
            for entry in ref["Value"]["StringWithMarkup"]:
                hphrase = entry["String"][0:4]
                if hphrase == 'Not ':
                    continue
                elif hphrase == 'Repo':
                    references[ref_id].companies = int(entry["String"].split()[8])
                    hazard_codes = [None]
                else:
                    hazard_codes.append(entry["String"][0:4])
            references[ref_id].hazard_codes = hazard_codes
        else:
            pass
            #print('\n ---NO GHS NOR ECHA---\n')
            #print(ref)

        if ref_name == "ECHA C&L Notifications Summary":
            #print('\n \x1b[1;32;80mECHA HERE\x1b[1;32;0m \n')
            #print(ref)
            entry = ref["Value"]["StringWithMarkup"][0]["String"]
            entry = entry.split()
            references[ref_id].companies = int(entry[5])
            references[ref_id].notifications = int(entry[8])
        else:
            pass
            #print('\n ---NO GHS NOR ECHA---\n')
            #print(ref)

    return references  # ,response_json

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 11:53:11 2022

@author: mpalermo
"""

import copy
from dataclasses import dataclass
from pubchemtools.communication import HTTP_request, PUBCHEM_load_balancer
from rdkit import Chem
import time


@dataclass
class ChemicalIDs:
    InChIKey: str = None
    InChI: str = None
    IUPAC: str = None
    SMILES: str = None
    PUBCHEM_CID: int = None


@dataclass
class GHSReferences:
    reference_id: str = None
    source_name: str = None
    hazard_codes: list = None
    companies: int = None
    notifications: int = None


@dataclass
class UserProperties:
    pass

    @property
    def count(self):
        return len(self.__dict__)

    def __setattr__(self, name, value, class_type=None):
        if class_type is not Compound:
            if name in self.__dict__:
                super().__setattr__(name, value)
            else:
                print("To set a new property, use the add_property instance method")
                return
        else:
            super().__setattr__(name, value)


class Compound:
    excluded_properties = ["_fetch_pairs"]
    builtin_properties = [
        "chemical_ids",
        "_vendors_data",
        "_vendors",
        "molecule",
        "_user_properties",
        "compound_in_PUBCHEM",
        "ghs_references",
        "_fetch_pairs",
        "_linked_libraries",
    ]

    def __init__(
        self,
        InChIKey=None,
        InChI=None,
        SMILES=None,
        IUPAC=None,
        sdf=None,
        mol=None,
        rdkitmol=None,
        autofetch=True,
    ):
        # Instance attributes
        self._vendors_data = None
        self._vendors = None
        self.molecule = None
        self._user_properties = UserProperties()
        self.compound_in_PUBCHEM = None  # will be set either to True or False
        self.ghs_references = {}
        self._linked_libraries = []

        # DEV: should be in try/except
        if sdf is not None:
            suppl = Chem.SDMolSupplier(sdf)
            self.molecule = suppl[0]
        if mol is not None:
            self.molecule = Chem.MolFromMolFile(mol)
        if rdkitmol is not None:
            self.molecule = rdkitmol

        if any([InChIKey, InChI, IUPAC, SMILES, sdf, mol, rdkitmol]):
            self.chemical_ids = ChemicalIDs(
                InChI=InChI, InChIKey=InChIKey, IUPAC=IUPAC, SMILES=SMILES
            )
        else:
            raise ValueError(
                "Compound object initialization requires at least one argument"
            )

        self._fetch_pairs = {
            self._vendors_url: self._set_vendors_pubresponse,
            self._ghs_url: self._set_ghs_pubresponse,
        }

        if autofetch:
            self.fetch_data()

    def _chemical_ids_url(self):

        # DEV: should add a try/except block here
        if self.molecule:
            self.chemical_ids.InChIKey = Chem.MolToInchiKey(self.molecule)

        chosen = None
        post = None

        if self.chemical_ids.InChIKey is not None:
            chosen = {"InChIKey": self.chemical_ids.InChIKey}
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{self.chemical_ids.InChIKey}/property/InChI,IsomericSMILES,IUPACName/JSON"
        elif self.chemical_ids.InChI is not None:
            chosen = {"InChI": self.chemical_ids.InChI}
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/property/InChIKey,IsomericSMILES,IUPACName/JSON"
            post = f"inchi={self.chemical_ids.InChI}"
        elif self.chemical_ids.SMILES is not None:
            chosen = {"IsomericSMILES": self.chemical_ids.SMILES}
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{self.chemical_ids.SMILES}/property/InChIKey,InChI,IUPACName/JSON"
        elif self.chemical_ids.IUPAC is not None:
            chosen = {"IUPACName": self.chemical_ids.IUPAC}
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/iupacname/{self.chemical_ids.IUPAC}/property/InChIKey,InChI,SMILES/JSON"

        return url, post

    def _vendors_url(self):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{self.chemical_ids.PUBCHEM_CID}/JSON/"
        post = None
        return url, post

    def _ghs_url(self):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{self.chemical_ids.PUBCHEM_CID}/JSON/?heading=GHS+Classification"
        post = None
        return url, post

    def _set_chemical_ids_pubresponse(self, response):
        identifiers = self.chemical_ids.__dict__

        # relies on definition order of dict - may be unreliable depending
        # on python base library implementation
        for identifier, value in identifiers.items():
            if value is not None:
                if identifier == "SMILES":
                    identifier = "IsomericSMILES"  # Pubchem distinguishes between isomeric and canonical SMILES
                chosen = {identifier: value}

        if response is None:
            print(
                f"No Chemical IDs were found for {list(chosen.keys())[0]}: {list(chosen.values())[0]}"
            )
        else:
            data = response["PropertyTable"]["Properties"][0] | chosen

            self.chemical_ids.PUBCHEM_CID = data["CID"]
            self.chemical_ids.InChIKey = data["InChIKey"]
            self.chemical_ids.InChI = data["InChI"]
            self.chemical_ids.SMILES = data["IsomericSMILES"]
            self.chemical_ids.IUPAC = data["IUPACName"]

            # Create molecular molfile
            self.molecule = Chem.MolFromInchi(self.chemical_ids.InChI)

    def _set_vendors_pubresponse(self, response):
        self._vendors_data = response["SourceCategories"]["Categories"][0]["Sources"]

    def _set_ghs_pubresponse(self, response):
        # data from Pubchem
        ref_sources = response["Record"]["Reference"]
        ref_data = response["Record"]["Section"][0]["Section"][0]["Section"][0][
            "Information"
        ]

        for source in ref_sources:
            ref_id = source["ReferenceNumber"]
            source_name = source["SourceName"]
            self.ghs_references[ref_id] = GHSReferences(
                reference_id=ref_id, source_name=source_name
            )

        for ref in ref_data:
            ref_id = ref["ReferenceNumber"]
            ref_name = ref["Name"]

            if ref_name == "GHS Hazard Statements":
                # print('\n \x1b[1;34;80mGHS Hazard Statements\x1b[1;34;0m \n')
                hazard_codes = []
                for entry in ref["Value"]["StringWithMarkup"]:
                    hphrase = entry["String"][0:4]
                    if hphrase == "Not ":
                        continue
                    elif hphrase == "Repo":
                        self.ghs_references[ref_id].companies = int(
                            entry["String"].split()[8]
                        )
                        hazard_codes = [None]
                    else:
                        hazard_codes.append(entry["String"][0:4])
                self.ghs_references[ref_id].hazard_codes = hazard_codes
            else:
                pass

            if ref_name == "ECHA C&L Notifications Summary":
                # print('\n \x1b[1;32;80mECHA HERE\x1b[1;32;0m \n')
                # print(ref)
                entry = ref["Value"]["StringWithMarkup"][0]["String"]
                entry = entry.split()
                self.ghs_references[ref_id].companies = int(entry[5])
                self.ghs_references[ref_id].notifications = int(entry[8])

    def fetch_data(self):
        # Chemical IDs
        url, post = self._chemical_ids_url()
        response = HTTP_request(url, post)
        self._set_chemical_ids_pubresponse(response)

        url, post = self._vendors_url()
        response = HTTP_request(url, post)
        self._set_vendors_pubresponse(response)

        url, post = self._ghs_url()
        response = HTTP_request(url, post)
        self._set_ghs_pubresponse(response)

    @property
    def vendors(self):
        return {vendor_entry["SourceName"] for vendor_entry in self._vendors_data}

    @property
    def user_properties(self):
        return self._user_properties

    def save(self, filename, file_format="sdf"):
        # DEV: should add restrains to arguments for file_format
        writer = Chem.SDWriter(filename + "." + file_format)

        molecule = copy.copy(self.molecule)
        molecule.SetProp("InChI", self.chemical_ids.InChI)
        molecule.SetProp("InChIKey", self.chemical_ids.InChIKey)
        molecule.SetProp("SMILES", self.chemical_ids.SMILES)
        molecule.SetProp("IUPAC", self.chemical_ids.IUPAC)

        molecule.SetProp("Vendors", ",".join(self.vendors))

        writer.write(molecule)

    def add_property(self, property_name, property_value):
        if property_name not in self._user_properties.__dict__:
            try:
                self._user_properties.__setattr__(
                    property_name, property_value, type(self)
                )
                for library in self._linked_libraries:
                    library._register_property(property_name)
            except Exception as e:
                print(str(e))
        else:
            self._user_properties.__setattr__(property_name, property_value, type(self))

    def __setattr__(self, name, value):
        if name not in self.__dict__ and name not in self.builtin_properties:
            print(name, value)
            print("To set a new property, use the add_property instance method")
            return
        else:
            super().__setattr__(name, value)


class Library:
    def __init__(self, compounds_list=None):

        self._registered_properties = set()

        if compounds_list is None:
            self.compounds_list = []
        else:
            self.compounds_list = compounds_list

        for compound in self.compounds_list:
            compound._linked_libraries.append(self)

    def fetch_data(self, parallel=True):
        # for compound in self.compounds_list:
        #     # using function as dict key to pair url request and setter!
        #     for url_fetcher, setter in compound._fetch_pairs.items():
        #         url, post = url_fetcher()
        #         response = HTTP_request(url, post)
        #         setter(response)
        initial_time = time.time()

        url_requests = []

        print("Fetching chemical IDs")
        # first, we need a CID...
        cid_requests = []
        for compound in self.compounds_list:
            cid_requests.append(compound._chemical_ids_url)

        # lets grab the chemical IDs from pubchem...
        data = PUBCHEM_load_balancer(cid_requests)

        # ... and set them in each compound object
        for compound in self.compounds_list:
            key = compound._chemical_ids_url
            try:
                compound._set_chemical_ids_pubresponse(data[key])
                compound.compound_in_PUBCHEM = True
            except:
                compound.compound_in_PUBCHEM = False

        print("\nFetching safety and vendors data")
        # Fetch other data
        for compound in self.compounds_list:
            for url_fetcher in compound._fetch_pairs.keys():
                if compound.compound_in_PUBCHEM:
                    url_requests.append(url_fetcher)

        data = PUBCHEM_load_balancer(url_requests)

        for compound in self.compounds_list:
            keys = compound._fetch_pairs.keys()
            for key in keys:
                setter = compound._fetch_pairs[key]
                try:
                    response = data[key]
                    setter(response)
                except:
                    pass
        print(f"Elapsed: ", time.time() - initial_time)

    def add(self, compounds_list):
        if type(compounds_list) is not list:
            raise TypeError(
                f"add method requires a list as input, while {type(compounds_list)} was provided."
            )
            return

        self.compounds_list.extend(compounds_list)

        for compound in compounds_list:
            for user_property in compound._user_properties.__dict__:
                self._register_property(user_property)

    def read(self, sdf_filepath):
        if sdf_filepath is not None:
            suppl = Chem.SDMolSupplier(sdf_filepath)
            # suppl = Chem.MultithreadedSDMolSupplier(sdf_filepath)

            for mol in suppl:
                self.compounds_list.append(Compound(rdkitmol=mol, autofetch=False))

    def save(self, filename, file_format="sdf"):

        sdf_out = Chem.SDWriter(filename)

        # output_molecules = []

        for compound in self.compounds_list:
            # avoids having replicated data in memory...
            if compound.molecule is not None:
                molecule = copy.copy(compound.molecule)

            else:
                molecule = Chem.rdchem.Mol()

            molecule.SetProp("InChI", str(compound.chemical_ids.InChI))
            molecule.SetProp("InChIKey", str(compound.chemical_ids.InChIKey))
            molecule.SetProp("SMILES", str(compound.chemical_ids.SMILES))
            molecule.SetProp("IUPAC", str(compound.chemical_ids.IUPAC))

            sdf_out.write(molecule)

        sdf_out.close()

    def _register_property(self, user_property):
        self._registered_properties.add(user_property)

    @property
    def in_pubchem(self):
        for compound in self.compounds_list:
            if compound.compound_in_PUBCHEM:
                yield compound

    def __iter__(self):
        for compound in self.compounds_list:
            yield compound

    def __getitem__(self, key):
        return self.compounds_list[key]

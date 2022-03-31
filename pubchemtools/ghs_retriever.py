#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 14:27:26 2021
@author: mattia
"""

from pubchemtools import pugrest_tools as prt
import requests
import threading, queue
import time
import datetime
import numpy as np


class CompoundLibrary:
    def __init__(self, compound_list, style="CAS"):
        """
        Parameters
        ----------
        compound_list : TYPE
            A list of compounds, either in CAS, smiles InChI or InChIKey format.
        style : TYPE, optional
            DESCRIPTION. The default is "CAS".
        Returns
        -------
        None.
        """

        self.references = {}  # dictionary containing compound key: dict of references
        self.results_queue = queue.Queue()  # self explanatory
        self._timings = (
            queue.Queue()
        )  # used to calculate the number of requests to PUBCHEM

        self._waiting_time = 1.9  # regulate speed

        # Retrieve references for the given library
        self._populate_references(compound_list, style)

        # self.clean_echa() # remove all sources and just leave the most reliable ECHA oource

    def _compound_to_ghs(self, compound, style="CAS", session=None):
        """
        Parameters
        ----------
        compound : str
            Smile/CAS/InChI/InChIKey of the compound.
        style : str, optional
            Specify the compound format. The default is "CAS".
        session : requests.sessions.Session, optional
            Session object to speed up API requests. The default is None.
        Returns
        -------
        references : dict
            Dict containing dicts of References objects.
        """

        # For regulating the number of requests to PUBCHEM
        start_time = time.time()

        # Choose the correct function depending on the compound format
        if style == "InChI":
            cid = prt.inchi_to_cid(compound, session)
        if style == "CAS":
            cid = prt.cas_to_cid(compound, session)
        if style == "InChIKey":
            cid = prt.inchikey_to_cid(compound, session)

        # Fetch references from PUBCHEM
        references = prt.cid_to_ghs(cid[0])

        # Store request end time and elapsed time
        self._timings.put([time.time(), time.time() - start_time, start_time])

        # Add references to the results queue
        self.results_queue.put([compound, references])

    def _populate_references(self, compound_list, style):
        session = requests.Session()

        # Queue of input compounds
        compounds_queue = queue.Queue()
        for compound in compound_list:
            compounds_queue.put(compound)

        self._start_time = time.time()  # for stat info total elapsed time

        # For displaying the total number of compounds
        self.compound_num = len(compound_list)

        # Counter to join groups of 50 threads
        threads = 0

        # Initialize to True so that program can start
        within_pubchem_limits = True

        while compounds_queue.empty() is False:
            new_threads = []  # group threads to be joined
            for thread_num in range(10):  # create 10 workers per second

                # Create workers only if we are not exceeding PUBCHEM volume limitations
                if (compounds_queue.empty() is False) and (within_pubchem_limits):
                    compound = compounds_queue.get()
                    worker = threading.Thread(
                        target=self._compound_to_ghs, args=[compound, style, session]
                    )
                    new_threads.append(worker)
                    worker.start()

            print(
                "\nProgress: ",
                str(self.results_queue.qsize()) + "/" + str(self.compound_num) + "\n",
            )

            within_pubchem_limits = (
                self._check_timings()
            )  # limits speed based on running time
            time.sleep(self._waiting_time)

            for thread in new_threads:
                thread.join()

        # Transform results queue to dictionary
        while self.results_queue.empty() is False:
            element = self.results_queue.get()
            compound = element[0]
            references = element[1]
            self.references[compound] = references

        print(
            "Total elapsed time: ",
            datetime.timedelta(seconds=time.time() - self._start_time),
        )

    def clean_echa(self):
        # Just keep the most relevant ECHA record based on number of companies

        compound_keys = list(self.references.keys())  # list of compounds

        # Go through all compounds
        for compound in compound_keys:
            compound_dict = self.references[compound]
            references_keys = list(compound_dict.keys())  # list of compound keys

            # Go through each reference for each molecule
            max_count = 0  # highest count of number of companies
            last_ref = None  # refid of the ref with highest number of companies

            # Cycle over reference for each compound
            for ref in references_keys:
                companies = compound_dict[ref].companies
                print("Compound ", compound)
                print("Reference ", ref, "Companies ", companies)

                # If company count for the ref is lower than another ref
                # or None (for non ECHA entries), delete entry
                if (companies is None) or (companies <= max_count):
                    print("deleting", max_count)
                    del self.references[compound][ref]  # continue
                elif companies >= max_count:
                    if last_ref is None:
                        print(
                            "Ref ",
                            ref,
                            " is the first entry with a populated companies field. Setting last_ref and max count.",
                        )
                        last_ref = ref
                        max_count = companies
                    elif last_ref is not None:
                        print(
                            "Deleting reference ",
                            last_ref,
                            " because ref ",
                            ref,
                            "has a higher count of companies",
                        )
                        del self.references[compound][last_ref]
                        last_ref = ref
                        max_count = companies
                        print("new last_ref", last_ref, " and max_count ", max_count)

    def _check_timings(self):
        queue_size = self._timings.qsize()
        # Do not check timing if we haven't started yet...
        if queue_size == 0:
            return True

        # use the last 500 records to evaluate timings
        if queue_size > 500:
            [self._timings.get() for i in range(queue_size - 500)]

        # Convert queue to array for manipulation
        timings = np.array(list(self._timings.queue))

        if queue_size > 10:
            req_per_sec = 10 / (time.time() - timings[-10, 2])
            # print('req per sec: ', req_per_sec)

            if req_per_sec < 5:
                self._waiting_time *= 0.95
            elif req_per_sec > 5:
                self._waiting_time *= 1.05
            # print('waiting time: ', self._waiting_time)

        elapsed_time = time.time() - timings[0, 0]
        running_time = sum(timings[:, 1])
        ratio = running_time / elapsed_time
        print("running time per elapsed time: ", ratio)

        if ratio > 6:  # empirically higher on purpose, as it leads closer to 5
            print("Regulating elapsed over running time")
            return False
        else:
            return True

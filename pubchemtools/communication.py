#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:22:38 2022

@author: mpalermo
"""

import requests
from requests.exceptions import HTTPError
import threading, queue
import time

import os

session = requests.Session()

# %% HTTP Catch errors and reuse session in HTTP request
def HTTP_request(url, post=None, other=None, attempts=3, queue=None):
    # if session is None:
    #     session = requests.Session()

    for i in range(1, attempts + 1):  # try three times

        try:
            response = session.post(url, data=post, timeout=2)
            # If the response was successful, no Exception will be raised
            #print(response.headers)
            response.raise_for_status()

        # DEV: should be completed with the full error list from PUBCHEM
        except HTTPError as http_err:
            if "PUGVIEW.NotFound" in str(http_err):
                # print(f"No GHS data found for compound CID {other}")
                return None

            elif "PUGREST.NotFound" in str(http_err):
                # print("PUGREST.NotFound")
                return None
            elif "PUGREST.ServerBusy" in str(http_err):
                print(f"\n PUBCHEM HTTP 503 error \"PUGREST.ServerBusy\" for url:\n  {url} \n Retrying: {i}/{attempts} attempts\n")
                time.sleep(1)
            else:
                print(str(http_err))
            break
            # print(f"HTTP error occurred: {http_err}")  # Python 3.6
        except Exception as err:
            print(f"Other error occurred: {err}")  # Python 3.6
            #print(f"\n Connection error. Retrying {other} - Attempt {i}/3")
            time.sleep(1)
        else:
            if queue is not None:
                queue.put([other, response.json()])
                return
            else:
                return response.json()
    else:
        print(f"Three attempts have failed for {url}, {post}")


def PUBCHEM_load_balancer(url_requests_list, total_urls=0, attempts=3):

    timings = [time.time()]
    initial_size = len(url_requests_list)
    last_size = 0

    def within_pubchem_limits(url_requests_list, current_time):
        nonlocal last_size
        size = initial_size - len(url_requests_list)
        if size % 1 == 0 and size != last_size:
            print(f"Completed {size}/{initial_size}")

        if size % 5 == 0 and size >= 5 and size != last_size:
            elapsed = current_time - timings[-1]

            # print(f"Last requests: {last_size}, current {size}")

            last_size = size
            timings.append(current_time)
            # print(f"Elapsed time: {elapsed}")

            speed = 5 / elapsed
            # start_time = time.time()

            if speed <= 5:
                return True
            else:
                print(f"Speed higher than 5 req/s: {speed}. Throttling.")
                time.sleep(1)
                return False
        return True

    url_requests_queue = queue.Queue()
    results_queue = queue.Queue()

    for url_request in url_requests_list:
        url_requests_queue.put(url_request)

    while url_requests_list:
        new_threads = []  # group threads to be joined
        for thread_num in range(5):  # create 10 workers per second
            if within_pubchem_limits(url_requests_list, time.time()):
                # try:
                #     url_request = url_requests_queue.get_nowait()
                # except queue.Empty:
                #     print("I have processed all urls")
                try:
                    url_request = url_requests_list.pop()
                except:
                    break

                url, post = url_request()
                # print("url: ", url)
                worker = threading.Thread(
                    target=HTTP_request, args=[url, post, url_request, attempts, results_queue]
                )
                new_threads.append(worker)
                worker.start()

        time.sleep(1.05)

        for thread in new_threads:
            thread.join()
            # print("joined")

    # print("Results queue size: ", results_queue.qsize())

    results = {}

    while not results_queue.empty():
        try:
            result = results_queue.get_nowait()
            results[result[0]] = result[1]

        except queue.Empty:
            pass

    return results

#!/usr/bin/env python3
"""fetch_pcnet2.py — fetches PCNet2.0 from NDEx and removes self-loop edges.
Source: 04_gwas.../1_NetColoc/1.1_OCD_dogCD_NetColoc_analysis_260521.ipynb

Reproduces the notebook's own interactome-fetch cell, not NetColoc's own netprop_zscore()
wrapper function: the notebook fetches the interactome and calls
netprop_zscore.calculate_heat_zscores() directly, never the wrapper, so the wrapper's own
internal cleanup step ("remove 'None' node") is not what runs here. Cleanup here is "remove
self-loop edges" (`G_int.remove_edges_from(nx.selfloop_edges(G_int))`), matching the notebook.

interactome_uuid defaults to PCNet2.0 (d73d6357-e87b-11ee-9621-005056ae23aa, 19,267 nodes /
3,852,119 edges).

The raw CX export for this network is >1GB — much larger than the summarized "19,267 nodes"
figure suggests, since CX format repeats node/edge attributes verbosely. The streamed download
can intermittently fail partway through (`IncompleteRead`/`ChunkedEncodingError`) for a transfer
this large — retried here with backoff rather than relying on Nextflow's own task-level retry,
which would re-run the entire task just to retry one HTTP call.
"""

import argparse
import pickle
import time

import networkx as nx
import ndex2
import requests

parser = argparse.ArgumentParser(description="Fetch PCNet2.0 (or another NDEx interactome) and remove self-loops.")
parser.add_argument("--uuid", default="d73d6357-e87b-11ee-9621-005056ae23aa", help="NDEx interactome UUID (default: PCNet2.0)")
parser.add_argument("--server", default="public.ndexbio.org", help="NDEx server")
parser.add_argument("--max-attempts", type=int, default=5, help="Retry attempts for the (large, occasionally flaky) network fetch")
parser.add_argument("out_file", help="Output pickle file for the fetched networkx.Graph")
args = parser.parse_args()

last_error = None
nice_cx = None
for attempt in range(1, args.max_attempts + 1):
    try:
        print(f"Fetch attempt {attempt}/{args.max_attempts}")
        nice_cx = ndex2.create_nice_cx_from_server(
            args.server,
            username=None,
            password=None,
            uuid=args.uuid,
        )
        break
    except (requests.exceptions.ChunkedEncodingError, requests.exceptions.ConnectionError) as e:
        last_error = e
        print(f"Attempt {attempt} failed: {e}")
        if attempt < args.max_attempts:
            time.sleep(10 * attempt)

if nice_cx is None:
    raise last_error

interactome = nice_cx.to_networkx()

interactome.remove_edges_from(nx.selfloop_edges(interactome))

print(f"Number of nodes: {len(interactome.nodes)}")
print(f"Number of edges: {len(interactome.edges)}")

with open(args.out_file, "wb") as f:
    pickle.dump(interactome, f)

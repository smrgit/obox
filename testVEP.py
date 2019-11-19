
import requests, sys
 
server = "https://rest.ensembl.org"
## ext = "/variant_recoder/human/AGT:c.803T>C?"
ext = "/variant_recoder/human/chr7:101837174:C:T"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
 r.raise_for_status()
 sys.exit()
 
decoded = r.json()
print(repr(decoded))



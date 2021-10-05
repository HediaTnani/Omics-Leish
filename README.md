# Omics-Leish

## Extracting metabolites, reactions and genes information from BIGG models

```
def extract_metabolites(x):
    a = []
    for m in x.metabolites:
        a.append(m.id)
    return(a)
```
```
def extract_reactions(x):
    b = []
    for r in x.reactions:
        b.append(r.id)
    return(b)
```
```
def extract_genes(x):
    c = []
    for g in x.genes:
        c.append(g.id)
    return(c)
```
## Parsing the Kegg API

```
def getkeggcompound(dono):
    from Bio.KEGG import REST
    import requests
    import json
    list_response = []
    result = []
    bigg2kegg = []
    l1 = []
    l2 = []
    for i in dono:  
        api_url =  "http://bigg.ucsd.edu/api/v2/universal/metabolites/"
        urls = api_url + i[:-2]
        #print(urls)
        response = requests.get(urls)
        if response.status_code == 200:
            list_response.append(response.json())
        else:
            pass
    for i in range(len(list_response)):
        #print(list_response[i])
        if 'bigg_id' in list_response[i]:
                #print(list_response[i]['bigg_id'])
                if 'KEGG Compound' in list_response[i]['database_links']:
                    l1.append(list_response[i]['bigg_id'])
                    l2.append(list_response[i]['database_links']['KEGG Compound'][0]['id'])
                    bigg2kegg.append(list_response[i]['bigg_id'] + ":" + list_response[i]['database_links']['KEGG Compound'][0]['id'])
                    result.append(REST.kegg_get("cpd:" + list_response[i]['database_links']['KEGG Compound'][0]['id'] ).read().split("\n"))


    d = dict(zip(l1, l2))
    conversion = pd.DataFrame.from_dict(d, orient="index").reset_index()
    conversion.columns = ['Bigg_id', 'Kegg_compound']
    return(conversion)
```


```
def parser(compound_file):
    import re
    pathways = []
    dictionary = {}
    current_section = None
    for line in compound_file.rstrip().split("\n"):
        #print(line)
        left_column = line[:12].strip()
        if left_column == "":
            pass
        else:
            current_section = left_column
        #print(current_section)
        if "PATHWAY" == current_section:
            right_column = line[12:]
            m = re.search("(.+?)\s+(.+)", right_column)
            #print(right_column)
            if m:
                dictionary[m.group(1)] = m.group(2)
    #print(dictionary)
    result = pd.DataFrame.from_dict(dictionary, orient="index").reset_index()
    result.columns = ['Map', 'Pathway']
    return result
```

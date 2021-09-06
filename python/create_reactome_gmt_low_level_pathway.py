import requests
from requests.api import get

reactomeGmtFile = "/home/spinicck/PhD/Data/gene-set/reactome_pathways_symbol.gmt"
newReactomeGmtFile = "/home/spinicck/PhD/Data/gene-set/reactome_pathways_only_lowest_level.gmt"

r = requests.get('https://reactome.org/ContentService/data/eventsHierarchy/9606')
humanEventsHierarchy = r.json()

def hasChildren (x):
    return ( 'children' in x.keys() ) and x['children']

def getOnlyLowestLevelPathwayId (event, ids):
    if (  type(event) is list ):
        for e in event:
            ids.update(getOnlyLowestLevelPathwayId(e, ids))
    elif (  type(event) is dict and hasChildren(event) ):
        hasNoPathway = True
        eventChildren = event.get('children')
        for e in eventChildren:
            if (e.get('type') == 'Pathway'):
                hasOnePathway = False
                ids.update(getOnlyLowestLevelPathwayId(e, ids))
        if (hasNoPathway):
            return [e.get('stId')]
    return ids
    


lowestLevelPathwayIds = getOnlyLowestLevelPathwayId (humanEventsHierarchy, set())

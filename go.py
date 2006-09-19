import _GOLib
import re
from sets import Set
from urllib import urlretrieve
import cPickle

data_dir=".//data"

#Currently loaded ontology (GeneOntologyDB)
loadedGO=None

#Currently loaded slim ontology (GeneOntologyDB)
loadedSlimGO=None

#Currently loaded annotation (AnnotationDB)
loadedAnnotation=None

namespaceDict={
    "biological_process":1, "P":1,
    "cellular_component":2, "C":2,
    "molecular_function":4, "F":4}

evidenceDict={"IMP":1, "IGI":2, "IPI":4, "ISS":8, "IDA":16, "IEP":32, "IEA":64,
              "TAS":128, "NAS":256, "ND":512, "IC":1024, "RCA":2048}

multiplicitySet=Set(["alt_id","is_a","subset","synonym","related_synonym","exact_synonym","broad_synonym","narrow_synonym",
                     "xref_analog","xref_unknown","relationship"])

annotationFields=["DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GOID", "DB_Reference","Evidence","With_From","Aspect",
                  "DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxon","Date","Assigned_by"]

def __evidenceToInt(evidenceList):
    if not evidenceList:
        return 4095
    evidence=0
    for e in evidenceList:
        if type(e)==str:
            e|=evidenceDict[e]
        elif type(e)==int:
            e|=e
    return evidence

def __evidenceToList(evidenceCode):
    evidence=[]
    for key, val in evidenceDict.items():
        if val&evidenceCode:
            evidence.append(key)
    return evidence

class AnnotationDB:
    annotation=None         #holds a C annotation structure do not touch
    geneNames=None          #a list of all gene names in annotation
    annotationList=None     #a list of all instances of class Annotation
    aliasMapper=None        #alias mapper maps known aliases to gene names (column3 DB_Object_Symbol of annotation file)
    geneNamesDict=None      #a dictionary mapping any known alias to a list [DB_Object_ID, DB_Object_Symbol [,DB_Object_Synonym]] i.d. all known names
    geneAnnotations=None    #a dictionary mapping gene name (DB_Object_Symbol) to a list of all instances of Annotation with this name

class GeneOntologyDB:
    ontology=None           #holds a C ontology structure do not touch
    terms=None              #holds a list of all instances of class GOTerm
    termDict=None           #a dictionary mapping GO id's to instances of class GOTerm
    termDescriptorDict=None #a dictionary mapping GO id's and alt_id's to a tuple (id, namespace, def, alt_id)
    aliasMapper=None        #a dictionary mapping alt_id's and id's to id's

class GOTerm:
    def __init__(self, term):
        self.original={}
        self.GOId=None
        self.parents=[]
        self.fullText=""
        self.aspect="unknown"
        self.alt_id=[]

        self.processBlock(term)

    def processBlock(self, text):
        for s in text.split("\n"):
            if ":" not in s:
                continue
            index=s.find(":")
            block=s[:index]
            body=s[index+1 :].strip(" ")
            if block in multiplicitySet:
                if self.__dict__.has_key(block):
                    self.__dict__[block].append(body)
                else:
                    self.__dict__[block]=[body]
            else:
                self.__dict__[block]=body            
        self.GOId=self.id
        self.parents=self.__dict__.get("is_a", [])
        self.parents=[s.split("!")[0].strip(" ") for s in self.parents]
        self.aspect=self.__dict__.get("namespace", "unknown")
        self.alt_id=self.__dict__.get("alt_id", [])
        self.fullText+=text

    def toTuple(self):
        return (self.GOId, self.parents)

class Annotation:
    def __init__(self, fullText):
        self.fullText=fullText
        self.original=fullText.split("\t")
        self.geneName=self.original[2]
        self.GOId=self.original[4].strip(" ")
        self.evidence=self.original[6].strip(" ")
        self.aspect=self.original[8].strip(" ")
        self.alias=self.original[10].split("|")
        for key, val in zip(annotationFields, self.original):
            self.__dict__[key]=val

    def toTuple(self):
        return (self.geneName, self.GOId, evidenceDict[self.evidence] , namespaceDict[self.aspect]) #evidenceMapper[self.evidence]

def setDataDir(datadir):
    """Set the data directory where the annotation and ontology files are stored (default:'.//data')"""
    global data_dir
    data_dir=datadir

def getDataDir():
    """Get the data directory where the annotation and ontology files are stored"""
    return data_dir

def loadAnnotation(organizm="sgd"):
    """Loads the annotation for the specified organizm"""
    global loadedAnnotation
    loadedAnnotation=loadAnnotationFrom(data_dir+"//gene_association."+organizm)#+".PyAnnotationDB")

def loadGO():
    """Loads the ontology from 'data_dir//gene_ontlogy.obo' where data_dir is set using setDataDir (default: .//data)"""
    global loadedGO
    loadedGO=loadOntologyFrom(data_dir+"//gene_ontology.obo")#.PyOntologyDB")
    
def mapTermId(TermId):
    """Maps the TermId to id if TermId a known alias for id (TermId can map to it self), if TermId unknown return None""" 
    return loadedGO.aliasMapper.get(TermId, None)

def mapGeneName(GeneName):
    """Maps the GeneName to name if GeneName a known alias for name (GeneName can map to it self), if GeneName unknown return None""" 
    return loadedAnnotation.aliasMapper.get(GeneName, None)

def GOTermFinder(clusterGeneList, referenceGenes=None, evidenceCodes=None, slimsOnly=False, aspect="P", progressCallback=None):
    """The method accepts a list of cluster genes, optionally a list of reference genes (otherwise all annotated genes appearing in the loaded annotation file are used),
    and optionally a list of annotation evidence codes to use, otherwise all evidence codes are used. The slimsOnly argument indicates if only GO slims are to be used,
    otherwise all GO terms are considered. The method returns a dictionary of significant GO terms, items are (cluster genes that map to each selected term, p-value,
    number of reference genes that map to each significant GO term). The progressCallback if present will be called with a single argument for all values in range [0..99]
    """
    if not referenceGenes:
        referenceGenes=loadedAnnotation.geneNames
    if slimsOnly and not loadedSlimGO:
        setSlims()
    annotation=loadedAnnotation.annotation
    goslim=loadedSlimGO and loadedSlimGO.ontology
    evidence=__evidenceToInt(evidenceCodes)
    aspect=namespaceDict[aspect]
    return _GOLib.GOTermFinder(clusterGeneList, referenceGenes, slimsOnly, evidence, aspect, annotation, loadedGO.ontology, goslim, progressCallback)

def findTerms(geneList, slimsOnly=False, directAnnotationOnly=False, evidenceCodes=None, reportEvidence=True):
    """For each gene in geneList search for matching GO terms. Argument slimsOnly restricts the GO terms to the slim set. The method returns a dictionary where key is a
    matching GO term and items are (gene, evidence) if reportEvidence == True [gene only, if reportEvidence=False] that map to the term. Climb the GO if directAnnotationOnly=False,
    otherwise report direct annotation only.
    """
    evidence=__evidenceToInt(evidenceCodes)
    if slimsOnly and not loadedSlimGO:
        setSlims()
        
    goslim=loadedSlimGO and loadedSlimGO.ontology
    annotation=loadedAnnotation.annotation
    result=_GOLib.findTerms(geneList, slimsOnly, directAnnotationOnly, evidence, reportEvidence, annotation, loadedGO.ontology, goslim)
        
    if reportEvidence:
        result=dict([(key, [(gene, __evidenceToList(evidence)) for gene ,evidence in val]) for key, val in result.items()])
    return result

def findGenes(GOTerms=[], directAnnotationOnly=False, evidenceCodes=None, reportEvidence=True, progressCallback=None):
    """Return a dictionary where key is a matching gene and items are (GO terms) or (GO term, list of evidences) from the GOterms list.
    (Note this will take a lot of time if the directAnnotationOnly=False)"""
    evidence=__evidenceToInt(evidenceCodes)
    result=_GOLib.findGenes(GOTerms, evidence, reportEvidence, directAnnotationOnly, loadedAnnotation.annotation, loadedGO.ontology, progressCallback)
    if reportEvidence:
        result=dict([(key,[(term, __evidenceToList(evidence)) for term, evidence in val]) for key, val in result.items()])
    return result        
    
def extractGODAG(GOTerms=[]):
    """Return the part of GO DAG that includes the listed GOTerms."""
    expanded=[]
    queue=GOTerms
    while queue:
        id=queue.pop()
        term=loadedGO.termDict.get(id, None)
        term=term or loadedGO.termDict.get(loadedGO.aliasMapper.get(id, None), None)
        if term and (term.id not in expanded):
            expanded.append(term.id)
            queue.extend(term.parents)
    terms=[loadedGO.termDict[id] for id in expanded if id in loadedGO.termDict]
    return terms      
            
def setSlims(slims=None):
    """Set the slim subset of a loaded ontology. slims can be:
        -a string identifier of a subsetdef: (e.g. "goslim_generic", "goslim_plant" ...)
        -a filename of a slim ontology (e.g. "goslim_generic.obo" ...)
        -a list of GO term id's
    """
    global loadedSlimGO
    goslims=[]
    slims= slims or "goslim_generic"
    try:
        loadedSlimGO=loadOntologyFrom(slims)
    except:
        if type(slims)==list:
            goslims=[loadedGO.termDict[term] for term in slims]
        else:
            goslims=filter(lambda t:slims in t.__dict__.get("subset",[]) , loadedGO.terms)
        loadedSlimGO=GeneOntologyDB()
        loadedSlimGO.ontology=_GOLib.parseGOTerms([g.toTuple() for g in goslims])
        loadedSlimGO.ontology.aliasMapper=loadedGO.aliasMapper
        loadedSlimGO.terms=goslims
        loadedSlimGO.termDict=loadedGO.termDict
        loadedSlimGO.aliasMapper=loadedGO.aliasMapper
        loadedSlimGO.termDescriptorDict=loadedGO.termDescriptorDict

def getSlims():
    """Returns the curently loaded slim terms"""
    return loadedSlimGO.terms

def downloadGO():
    """Downloads the curent gene ontology from http://www.geneontology.org/ontology/gene_ontology.obo"""
    urlretrieve("http://www.geneontology.org/ontology/gene_ontology.obo", data_dir+"//gene_ontology.obo")
    file=open(data_dir+"//gene_ontology.obo")
    data=file.read()
    c=re.compile("\[Term\].*?\n\n",re.DOTALL)
    match=c.findall(lines)
    go=parseGeneOntology(match)
    cPickle.dump(go, open(data_dir+"gene_ontology.obo.PyOntologyDB", "w"))

def downloadAnnotation(organizm="sgd"):
    """Downloads the annotation for the specified organizm (e.g. "sgd", "fb", "mgi",...)"""
    urlretrieve("http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association."+organizm+".gz",
                data_dir+"//gene_association."+organizm+".gz")
    from gzip import GzipFile
    gfile=GzipFile(data_dir+"//gene_association."+organizm+".gz","r")
    data=gfile.readlines()
    file=open(data_dir+"//gene_association."+organizm,"w")
    file.writelines(data)
    __splitAnnotation(data, organizm)
    anno=parseAnnotation(data)
    import cPickle
    cPickle.dump(anno, open(data_dir+"//gene_association."+organizm+".PyAnnotationDB", "w"))

def listOrganizms():
    """Connect to http://www.geneontology.org/GO.current.annotations.shtml, parse out the organism names
    appearing in the table, and return the list of organisms."""
    try:
        urlretrieve("http://www.geneontology.org/GO.current.annotations.shtml", data_dir+"//annotations.shtml")
    except:
        print "Failed to connect to http://www.geneontology.org/GO.current.annotations.shtml. Trying to find a local copy"
    file=open(data_dir+"//annotations.shtml")
    data=file.read()
    match=re.findall(r'http://www\.geneontology\.org/cgi-bin/downloadGOGA\.pl/gene_association\..+?gz', data)

    organizms=[]
    for m in match:
        #names=re.findall(r'>.+?<', m)
        organizms.append((m.split(".")[-2]))
    return organizms

def parseGeneOntology(data):
    terms=[]
    termDict={}
    aliasMapper={}
    goTermDict={}
    for term in data:
        t=GOTerm(term)
        termDict[t.id]=t
        terms.append(t)
        alt=t.__dict__.get("alt_id", [])
        aliasMapper.update(dict([(alt_id.strip(" "), t.id) for alt_id in alt]))
        aliasMapper[t.id]=t.id
        tt=(t.id, t.namespace, t.__dict__.get("def",""), alt)
        for alt_id in alt:
            goTermDict[alt_id]=tt
        goTermDict[t.id]=tt
    GO=GeneOntologyDB()
    GO.terms=terms
    GO.termDict=termDict
    GO.aliasMapper=aliasMapper
    GO.termDescriptorDict=goTermDict
    return GO

def loadOntologyFrom(filename):
    if filename.endswith(".PyOntologyDB"):
        db=cPickle.load(open(filename))
    else:
        file=open(filename)
        data=file.read()
        c=re.compile("\[Term\].*?\n\n",re.DOTALL)
        match=c.findall(data)
        db=parseGeneOntology(match)
    db.ontology=_GOLib.parseGOTerms([t.toTuple() for t in db.terms])
    db.ontology.aliasMapper=db.aliasMapper
    db.__file__=filename
    return db

def parseAnnotation(data):
    aliasMapper={}
    annotationList=[]
    geneNames=Set()
    geneNamesDict={}
    geneAnnotation={}
    for line in data:
        if line.startswith("!"):
            continue
        a=Annotation(line)
        if not a.geneName or not a.GOId:
            continue
        if a.geneName not in geneNames:
            geneNames.add(a.geneName)
            geneAnnotation[a.geneName]=[a]
            for alias in a.alias:
                aliasMapper[alias]=a.geneName
            aliasMapper[a.geneName]=a.geneName
            names=[a.original[1], a.original[2]]
            names.extend(a.alias)
            for n in names:
                geneNamesDict[n]=names
        else:
            geneAnnotation[a.geneName].append(a)
        annotationList.append(a)
    a=AnnotationDB()
    a.annotationList=annotationList
    a.aliasMapper=aliasMapper
    a.geneNames=list(geneNames)
    a.geneNamesDict=geneNamesDict
    a.geneAnnotations=geneAnnotation
    return a

def loadAnnotationFrom(filename):
    if filename.endswith(".PyAnnotationDB"):
        anno=cPickle.load(open(filename))
    else:
        file=open(filename)
        data=file.readlines()
        anno=parseAnnotation(data)
        
    anno.annotation=_GOLib.parseAnnotation([a.toTuple() for a in anno.annotationList if "NOT" not in a.Qualifier])
    anno.annotation.aliasMapper=anno.aliasMapper
    anno.__file__=filename
    return anno
    
def __test():
    def call(i): print i
    setDataDir("E://GO//data")
    print "Loading GO"
    loadGO()
    print "Loading annotation"
    loadAnnotation()
    print "Loading slims"
    setSlims()
    print "Finding terms"
    print GOTermFinder(loadedAnnotation.geneNames[:20], progressCallback=call)
    print "Finding slim terms"
    print GOTermFinder(loadedAnnotation.geneNames[20:30], slimsOnly=True, aspect="F")
    print "Finding terms"
    print findTerms(loadedAnnotation.geneNames[100:120])
    print "Finding direct slim terms"
    terms=findTerms(loadedAnnotation.geneNames[200:210], slimsOnly=True, directAnnotationOnly=True)
    print terms
    print "Extracting GO Dag"
    print extractGODAG(terms.keys())
    print "Finding genes"
    print findGenes(terms.keys()[:min(len(terms.keys()),3)], progressCallback=call)
    print findGenes(["GO:0005763","GO:0003735","GO:0042255","GO:0043037"], progressCallback=call)
    
if __name__=="__main__":
    __test()
    
    
"""A library for handling gene ontologies.
"""
import _GOLib
import re
from sets import Set
from urllib import urlretrieve
import cPickle
import os

prefix=os.path.split(__file__)[0]
if prefix:
    data_dir=prefix+"//data"
else:
    data_dir=".//data"
    
#Currently loaded ontology (GeneOntologyDB)
loadedGO=None

#Currently loaded slim ontology (GeneOntologyDB)
loadedSlimGO=None

#Currently loaded annotation (AnnotationDB)
loadedAnnotation=None

#A dictionary for mapping gene aliases
geneMapper={}

#A dictionary for mapping GOIs's
termMapper={}

namespaceDict={
    "biological_process":1, "P":1,
    "cellular_component":2, "C":2,
    "molecular_function":4, "F":4}

evidenceDict={"IMP":1, "IGI":2, "IPI":4, "ISS":8, "IDA":16, "IEP":32, "IEA":64,
              "TAS":128, "NAS":256, "ND":512, "IC":1024, "RCA":2048}

evidenceTypesOrdered = [
'IMP',
'IGI',
'IPI',
'ISS',
'IDA',
'IEP',
'IEA',
'TAS',
'NAS',
'ND',
'IC'
]

evidenceTypes = {
'IMP': 'inferred from mutant phenotype',
'IGI': 'inferred from genetic interaction', ## [with <database:gene_symbol[allele_symbol]>]',
'IPI': 'inferred from physical interaction', ## [with <database:protein_name>]',
'ISS': 'inferred from sequence similarity', ## [with <database:sequence_id>] ',
'IDA': 'inferred from direct assay',
'IEP': 'inferred from expression pattern',
'IEA': 'inferred from electronic annotation', ## [to <database:id>]',
'TAS': 'traceable author statement',
'NAS': 'non-traceable author statement',
'ND': 'no biological data available ',
'IC': 'inferred by curator'
}

multiplicitySet=Set(["alt_id","is_a","subset","synonym","related_synonym","exact_synonym","broad_synonym","narrow_synonym",
                     "xref_analog","xref_unknown","relationship"])

annotationFields=["DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GOID", "DB_Reference","Evidence","With_From","Aspect",
                  "DB_Object_Name","DB_Object_Synonym","DB_Object_Type","taxon","Date","Assigned_by"]

annotationFieldsDict={"DB":0,
                      "DB_Object_ID":1,
                      "DB_Object_Symbol":2,
                      "Qualifier":3,
                      "GO_ID":4,
                      "GO ID":4,
                      "DB_Reference":5,
                      "DB:Reference":5,
                      "Evidence_code":6,
                      "Evidence code":6,
                      "With_or_From":7,
                      "With (or) From":7,
                      "Aspect":8,
                      "DB_Object_Name":9,
                      "DB_Object_Synonym":10,
                      "DB_Object_Type":11,
                      "taxon":12,
                      "Date":13,
                      "Assigned_by":14}

def __evidenceToInt(evidenceList):
    if not evidenceList:
        return 4095
    evidence=0
    for e in evidenceList:
        if type(e)==str:
            evidence|=evidenceDict[e]
        elif type(e)==int:
            evidence|=e
    return evidence

def __evidenceToList(evidenceCode):
    evidence=[]
    for key, val in evidenceDict.items():
        if val&evidenceCode:
            evidence.append(key)
    return evidence

class AnnotationDB:
    """An object holding the annotation database.
    members:
        -geneNames      -- Names of all the genes in the annotation
        -annotationList -- List of instances of Annotation class holding the details for each annotation record
        -aliasMapper    -- Alias mapper maps known aliases to gene names (column3 DB_Object_Symbol of annotation file)
        -geneNamesDict  -- A dictionary mapping any known alias to a list [DB_Object_ID, DB_Object_Symbol [,DB_Object_Synonym]] i.d. all known names
        -geneAnnotations-- A dictionary mapping gene name (DB_Object_Symbol) to a list of all instances of Annotation with this name
    """
    __annotation=None       #holds a C annotation structure for internal use(do not touch!)
    geneNames=None          #a list of all gene names in annotation
    annotationList=None     #a list of all instances of class Annotation
    aliasMapper=None        #alias mapper maps known aliases to gene names (column3 DB_Object_Symbol of annotation file)
    geneNamesDict=None      #a dictionary mapping any known alias to a list [DB_Object_ID, DB_Object_Symbol [,DB_Object_Synonym]] i.d. all known names
    geneAnnotations=None    #a dictionary mapping gene name (DB_Object_Symbol) to a list of all instances of Annotation with this name

class GeneOntologyDB:
    """An object holding the ontology database.
    members:
        -terms              - List of instances of GOTerm class holding the details for each term
        -termDict           - A dictionary mapping GO id's to an instance of GOTerm class with that id
        -termDescriptorDict - A dictionary mapping GO id's and alt_id's to a tuple (id, namespace, def, alt_id)
        -aliasMapper        - A dictionary mapping alt_id's and id's to id's
    """
    __ontology=None         #holds a C ontology structure for internal use(do not touch!)
    terms=None              #holds a list of all instances of class GOTerm
    termDict=None           #a dictionary mapping GO id's to instances of class GOTerm
    termDescriptorDict=None #a dictionary mapping GO id's and alt_id's to a tuple (id, namespace, def, alt_id)
    aliasMapper=None        #a dictionary mapping alt_id's and id's to id's

class GOTerm:
    """Holds the data for one term read from the ontology file. All fields from the ontology can be accsessed by their
    original name (e.g. the is_a field in the ontology can be accsessed like term.is_a), except for the def field that
    interferes with the python def statment and can be accesed like term._def or term.__dict__['def'].
    The fields that can have multiple values are stored as lists of strings otherwise the string after the field name from the
    ontology is supplied. If a field is not defined accessing it will result in an exception.
    The object also provides the folowing data memebers for quicker accsess: GOId, aspect, parents(list of GOId's of parents terms)."""
    def __init__(self, term):
        self.original={}
        self.GOId=None
        self.parents=[]
        self.rType={} #relationship type with a parent
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
        for p in self.parents:
            self.rType[p]="is_a"
        for r in self.__dict__.get("relationship",[]):
            s=r.split("!")[0].strip(" ").split(" ")
            rT,id=s[0],s[1]
            self.parents.append(id)
            self.rType[id]=rT
        self.aspect=self.__dict__.get("namespace", "unknown")
        self.alt_id=self.__dict__.get("alt_id", [])
        self.fullText+=text

    def __getattr__(self, name):
        if name=="_def":
            return self.__dict__["def"]
        else:
            raise AttributeError(name)
        
    def toTuple(self):
        return (self.GOId, self.parents)

class Annotation:
    """Holds the data for an annotation record read from the annotation file. Fields can be
    accessed with the names: DB, DB_Object_ID, DB_Object_Symbol, Qualifier, GO_ID, DB_Reference,
    Evidence_code, With_or_From, Aspect, DB_Object_Name, DB_Object_Synonym, DB_Object_Type, taxon,
    Date, Assigned_by (e.g. rec.GO_ID)
    or by supplying the original name of the field (see http://geneontology.org/GO.annotation.shtml#file)
    to the get method (e.g. rec.get("GO ID"))
    The object also provides the folowing data members for quicker access: geneName, GOId, evidence,
    aspect and alias(a list of aliases)
    """
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

    def __getattr__(self, name):
        if name in annotationFieldsDict:
            return self.original[annotationFieldsDict[name]]
        else:
            raise AttributeError(name)

    def get(self, name):
        if name in annotationFieldsDict:
            return self.original[annotationFieldsDict[name]]
        else:
            raise ValueError(name)
        
    def toTuple(self):
        return (self.geneName, self.GOId, evidenceDict[self.evidence] , namespaceDict[self.aspect]) #evidenceMapper[self.evidence]

def setDataDir(datadir):
    """Set the data directory where the annotation and ontology files are stored (by default a directory named data
    located where the GOLib is instaled e.g. ...site-packages/GOLib/data)"""
    global data_dir
    data_dir=datadir

def getDataDir():
    """Get the data directory where the annotation and ontology files are stored (by default a directory named data
    located where the GOLib is instaled e.g. ...site-packages/GOLib/data)"""
    return data_dir

def loadAnnotation(organism="sgd", forceReload=False, progressCallback=None):
    """Loads the annotation for the specified organism"""
    global loadedAnnotation
    if loadedAnnotation and loadedAnnotation.__file__.endswith(organism) and not forceReload:
        return
    loadedAnnotation=loadAnnotationFrom(data_dir+"//gene_association."+organism, progressCallback)#+".PyAnnotationDB")
    global geneMapper
    geneMapper=loadedAnnotation.aliasMapper

def loadGO(forceReload=False, progressCallback=None):
    """Loads the ontology from 'data_dir//gene_ontlogy.obo' where data_dir is set using setDataDir (default: .//data)"""
    global loadedGO
    if loadedGO and not forceReload:
        return
    loadedGO=loadOntologyFrom(data_dir+"//gene_ontology.obo", progressCallback)#.PyOntologyDB")
    global termMapper
    termMapper=loadedGO.aliasMapper
    
def mapTermId(TermId):
    """Maps the TermId to id if TermId a known alias for id (TermId can map to it self), if TermId unknown return None""" 
    return loadedGO.aliasMapper.get(TermId, None)

def mapGeneName(GeneName):
    """Maps the GeneName to name if GeneName a known alias for name (GeneName can map to it self), if GeneName unknown return None""" 
    return loadedAnnotation.aliasMapper.get(GeneName, None)

def mapGeneNames(names=[]):
    return filter(lambda a:a, map(mapGeneName, names))

def __filterSlimsGOId(d):
    slims=[t.GOId for t in getSlims()]
    slims=filter(lambda id:id in slims, d.keys())
    return dict([(id, d[id]) for id in slims])

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
    #clusterGeneList=mapGeneNames(clusterGeneList)
    #referenceGenes=mapGeneNames(referenceGenes)
        
    annotation=loadedAnnotation.__annotation
    goslim=loadedSlimGO and loadedSlimGO.__ontology
    evidence=__evidenceToInt(evidenceCodes)
    aspect=namespaceDict[aspect]
    result=_GOLib.GOTermFinder(clusterGeneList, referenceGenes, slimsOnly, evidence, aspect, annotation, loadedGO.__ontology, goslim, progressCallback)
    if slimsOnly:
        return __filterSlimsGOId(result)
    else:
        return result
    
def findTerms(geneList, slimsOnly=False, aspect=["F","C","P"], directAnnotationOnly=False, evidenceCodes=None, reportEvidence=True, progressCallback=None):
    """For each gene in geneList search for matching GO terms. Argument slimsOnly restricts the GO terms to the slim set. The method returns a dictionary where key is a
    matching GO term and items are (gene, evidence) if reportEvidence == True [gene only, if reportEvidence=False] that map to the term. Climb the GO if directAnnotationOnly=False,
    otherwise report direct annotation only.
    """
    evidence=__evidenceToInt(evidenceCodes)
    if slimsOnly and not loadedSlimGO:
        setSlims()
    aa=0
    for a in aspect:
        aa|=namespaceDict[a]
    #geneList=mapGeneNames(geneList)
    goslim=loadedSlimGO and loadedSlimGO.__ontology
    annotation=loadedAnnotation.__annotation
    result=_GOLib.findTerms(geneList, slimsOnly, directAnnotationOnly, aa, evidence, reportEvidence, annotation, loadedGO.__ontology, goslim, progressCallback)
    if slimsOnly:
        result=__filterSlimsGOId(result)
    if reportEvidence:
        result=dict([(key, [(gene, __evidenceToList(evidence)) for gene ,evidence in val]) for key, val in result.items()])
    return result

def findGenes(GOTerms=[], directAnnotationOnly=False, evidenceCodes=None, reportEvidence=True, progressCallback=None):
    """Return a dictionary where key is a matching gene and items are (GO terms) or (GO term, list of evidences) from the GOterms list.
    (Note this will take a lot of time if the directAnnotationOnly=False)"""
    evidence=__evidenceToInt(evidenceCodes)
    result=_GOLib.findGenes(GOTerms, evidence, reportEvidence, directAnnotationOnly, loadedAnnotation.__annotation, loadedGO.__ontology, progressCallback)
    if reportEvidence:
        result=dict([(key,[(term, __evidenceToList(evidence)) for term, evidence in val]) for key, val in result.items()])
    return result        
    
def extractGODAG(GOTerms=[]):
    """Return the part of GO DAG that includes the listed GOTerms."""
    expanded=[]
    queue=list(GOTerms)
    while queue:
        id=queue.pop()
        term=loadedGO.termDict.get(id, None)
        term=term or loadedGO.termDict.get(loadedGO.aliasMapper.get(id, None), None)
        if term and (term.id not in expanded):
            expanded.append(term.id)
            queue.extend(term.parents)
    terms=[loadedGO.termDict[id] for id in expanded if id in loadedGO.termDict]
    return terms      

def __DAGDepth(term, cache={}):
	if term.parents:
		d=[]
		for t in term.parents:
			if t in cache:
				d.append(cache[t]+1)
			else:
				d.append(__DAGDepth(loadedGO.termDict[t], cache)+1)
		depth=min(d)
		cache[term.id]=depth
		return depth
	else:
		return 1

def DAGDepth(DAGTerms=[]):
    """returns the maximum depth of terms in DAGTerms"""
    cache={}
    #DAGTerms=[type(t)==GOTerm and t or loadedGO.termDict[t] for t in DAGTerms]
    return max([__DAGDepth(t, cache) for t in DAGTerms])

def DAGFilterForDepth(DAGTerms, depth):
    cache={}
    return filter(lambda t:__DAGDepth(t,cache)<=depth, DAGTerms)

def extractGODAGToDisplay(geneList):
    """returns the part of the GO that all of the genes in geneList map to
    """
    terms=findTerms(geneList)
    return terms.keys()
    
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
        loadedSlimGO.__ontology=_GOLib.parseGOTerms([g.toTuple() for g in goslims], loadedGO.aliasMapper)
        #loadedSlimGO.__ontology.aliasMapper=loadedGO.aliasMapper
        loadedSlimGO.terms=goslims
        loadedSlimGO.termDict=loadedGO.termDict
        loadedSlimGO.aliasMapper=loadedGO.aliasMapper
        loadedSlimGO.termDescriptorDict=loadedGO.termDescriptorDict

def getSlims():
    """Returns the curently loaded slim terms"""
    return loadedSlimGO.terms

class __progressCallWrapper:
        def __init__(self,callback):
            self.callback=callback
        def __call__(self, bCount, bSize, fSize):
            #print bCount, bSize, fSize
            if fSize==-1:
                fSize=10000000
            self.callback(100*bCount*bSize/fSize)
            
def downloadGO(progressCallback=None):
    """Downloads the curent gene ontology from http://www.geneontology.org/ontology/gene_ontology.obo"""
    urlretrieve("http://www.geneontology.org/ontology/gene_ontology.obo", data_dir+"//gene_ontology.obo", progressCallback and __progressCallWrapper(progressCallback))
    file=open(data_dir+"//gene_ontology.obo")
    data=file.read()
    c=re.compile("\[Term\].*?\n\n",re.DOTALL)
    match=c.findall(data)
    go=parseGeneOntology(match)
    #cPickle.dump(go, open(data_dir+"gene_ontology.obo.PyOntologyDB", "w"))
    
def downloadAnnotation(organism="sgd", progressCallback=None):
    """Downloads the annotation for the specified organism (e.g. "sgd", "fb", "mgi",...)"""
        
    #urlretrieve("http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association."+organism+".gz",
    #            data_dir+"//gene_association."+organism+".gz", progressCallback and __progressCallWrapper(progressCallback))
    urlretrieve("http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association."+organism+".gz?rev=HEAD",
                data_dir+"//gene_association."+organism+".gz", progressCallback and __progressCallWrapper(progressCallback))
    from gzip import GzipFile
    gfile=GzipFile(data_dir+"//gene_association."+organism+".gz","r")
    data=gfile.readlines()
    file=open(data_dir+"//gene_association."+organism,"w")
    file.writelines(data)
    #__splitAnnotation(data, organism)
    anno=parseAnnotation(data)
    import cPickle
    cPickle.dump(anno.aliasMapper.keys(), open(data_dir+"//gene_names."+organism, "w"))

def getCachedGeneNames(organism="sgd"):
    import cPickle
    return cPickle.load(open(data_dir+"//gene_names."+organism))

def listOrganisms():
    """Connect to http://www.geneontology.org/GO.current.annotations.shtml, parse out the organism names
    appearing in the table, and return the list of organisms."""
    try:
        urlretrieve("http://www.geneontology.org/GO.current.annotations.shtml", data_dir+"//annotations.shtml")
    except:
        print "Failed to connect to http://www.geneontology.org/GO.current.annotations.shtml. Trying to find a local copy"
    file=open(data_dir+"//annotations.shtml")
    data=file.read()
    #match=re.findall(r'http://www\.geneontology\.org/cgi-bin/downloadGOGA\.pl/gene_association\..+?gz', data)
    match=re.findall(r'http://cvsweb\.geneontology\.org/cgi-bin/cvsweb\.cgi/go/gene-associations/gene_association\..+?gz\?rev=HEAD', data)
    organisms=[]
    for m in match:
        #print m
        #names=re.findall(r'>.+?<', m)
        organisms.append(m.split(".")[-2])
    return organisms

def listDownloadedOrganisms():
    """Returns a list with organism names off all local downloaded annotations"""
    import os
    files=os.listdir(data_dir)
    files=filter(lambda n: n.startswith("gene_association") and not n.endswith(".gz"), files)
    return [s[17:] for s in files]

def parseGeneOntology(data, progressCallback=None):
    terms=[]
    termDict={}
    aliasMapper={}
    goTermDict={}
    datalen=len(data)
    milestones=Set(range(0,datalen,datalen/10))
    for i, term in enumerate(data):
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
        if progressCallback and i in milestones:
            progressCallback(100.0*i/datalen)
    GO=GeneOntologyDB()
    GO.terms=terms
    GO.termDict=termDict
    GO.aliasMapper=aliasMapper
    GO.termDescriptorDict=goTermDict
    return GO

def loadOntologyFrom(filename,progressCallback=None):
    if filename.endswith(".PyOntologyDB"):
        db=cPickle.load(open(filename))
    else:
        file=open(filename)
        data=file.read()
        c=re.compile("\[Term\].*?\n\n",re.DOTALL)
        match=c.findall(data)
        db=parseGeneOntology(match, progressCallback)
    db.__ontology=_GOLib.parseGOTerms([t.toTuple() for t in db.terms], db.aliasMapper)
    #db.__ontology.aliasMapper=db.aliasMapper
    db.__file__=filename
    return db

def parseAnnotation(data, progressCallback=None):
    aliasMapper={}
    annotationList=[]
    geneNames=Set()
    geneNamesDict={}
    geneAnnotation={}
    datalen=len(data)
    milestones=Set(range(0,datalen, datalen/10))
    for i,line in enumerate(data):
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
        if progressCallback and i in milestones:
            progressCallback(100.0*i/datalen)
    a=AnnotationDB()
    a.annotationList=annotationList
    a.aliasMapper=aliasMapper
    a.geneNames=list(geneNames)
    a.geneNamesDict=geneNamesDict
    a.geneAnnotations=geneAnnotation
    return a

def loadAnnotationFrom(filename, progressCallback=None):
    if filename.endswith(".PyAnnotationDB"):
        anno=cPickle.load(open(filename))
    else:
        file=open(filename)
        data=file.readlines()
        anno=parseAnnotation(data, progressCallback)
        
    anno.__annotation=_GOLib.parseAnnotation([a.toTuple() for a in anno.annotationList if "NOT" not in a.Qualifier], anno.aliasMapper)
    #anno.__annotation.aliasMapper=anno.aliasMapper
    anno.__file__=filename
    return anno
    
def __getattr__(self, name):
	if name=="geneMapper":
		return loadedAnnotation.aliasMapper
	elif name=="termMapper":
		return loadedGO.aliasMapper
	else:
		raise ValueError(name)

def __test():
    def call(i): print i
    setDataDir("E://orangecvs//GOLib//data")
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
    print findGenes(terms.keys()[:min(len(terms.keys()),3)], progressCallback=call)#,directAnnotationOnly=True)
    print findGenes(["GO:0005763","GO:0003735","GO:0042255","GO:0043037"], progressCallback=call)#,directAnnotationOnly=True)

if not listDownloadedOrganisms():
    print "Warning!!! No downloaded annotations found!!!"
    print "You can download annotations using the downloadAnnotation function."
    print "e.g. go.downloadAnnotation('sgd')"
try:
    open(data_dir+"//gene_ontology.obo")
except:
    print "Warning!!! No downloaded ontology found!!!"
    print "You can download it using the downloadGO function."
    print "e.g. go.downloadGO()"
    
    
if __name__=="__main__":
    __test()
    
    
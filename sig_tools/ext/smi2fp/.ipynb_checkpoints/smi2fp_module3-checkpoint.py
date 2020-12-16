"""
A  module with functions to query CMap Mongo DB,
convert SMILES strings to binary fingerprints using pybel,
and read and write JSON files.

Version 1.0
E-mail questions to: morzech@broadinstitute.org
Marek Orzechowski 08/19/2015

"""

import sys
import ssl
import copy
import urllib
import json
import io
import csv
from re import match
from openbabel import pybel
from rdkit.Chem import AllChem
# For drawing molecules
from rdkit.Chem import Draw
import numpy
import pandas as pd
#import cmap.io.gct as gct
import cmapPy.pandasGEXpress.GCToo as gctoo
import cmapPy.pandasGEXpress.write_gctx as wgx
import cmapPy.pandasGEXpress.write_gct as wg

### function to query MongoDB and to write a json file with results of the query
def fetch_json(url, outfile):
    # a workaround for the [SSL: CERTIFICATE_VERIFY_FAILED]
    if hasattr(ssl, '_create_unverified_context'):
        ssl._create_default_https_context = ssl._create_unverified_context

    urlstring = urllib.urlopen(url).read()
    urljson = json.loads(urlstring)
    write_json(urljson, 'trt_cp.json')
    return urljson


### function to write a file in a JSON format to a disk
def write_json(json_obj, outfile):
    with io.open(outfile, 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(json_obj, ensure_ascii=False)))


### function to read a JSON file from a disk
def read_json(infile):
    with io.open(infile, 'r', encoding='utf-8') as f:
        dataJSON = json.load(f, encoding='utf-8')
    return dataJSON


### function that checks for missing SMILES and removes such elements from the query result
def check_missing(query_list):
    # initialize counters of all elements of the list with present (p) and missing (m) keys
    p = 0;
    m = 0;
    elements_to_remove = []
    if type(query_list) == dict:
        pass
    elif type(query_list) == list:
        for i in range(len(query_list)):
            if 'canonical_smiles' in query_list[i]:
                p = p + 1
                if query_list[i]['canonical_smiles'] == "-666":
                    print('SMILES is -666 in {0:s}'.format(query_list[i]['pert_id']))
                    elements_to_remove.append(i)
                    p = p - 1
                    m = m + 1
            else:
                m = m + 1
                print('SMILES is missing in {0:s}'.format(query_list[i]['pert_id']))
                elements_to_remove.append(i)

    print("Present: {0:d}, Missing: {1:d}, Total: {2:d}".format(p, m, p + m))
    # print "Elements to remove: %s" % elements_to_remove
    # remove missing elements
    for i in list(reversed(elements_to_remove)):
        query_list.pop(i)

    print(len(query_list))


# return query_list;

### draw structure from a SMILES string
def smiles2image(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    # PIL image, save using img.save(filename)
    return img
    
### function to convert SMILES to binary fingerprints (uses openbabel)
def smiles2fpt(smiles, fptt='FP2', radius=2, bit_size=1024):
    try:
        if fptt == 'FCFP':
            mol = AllChem.MolFromSmiles(smiles)
            fps = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=bit_size, useFeatures=True)
            fpt = [int(x) for x in list(fps.ToBitString())]
        elif fptt == 'ECFP':
            mol = AllChem.MolFromSmiles(smiles)
            # Note ECFP4 is diameter=4, which corresponds to MorganFP radius=2
            fps = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=bit_size, useFeatures=False)
            fpt = [int(x) for x in list(fps.ToBitString())]
        elif fptt == 'RDKFP':
            mol = AllChem.MolFromSmiles(smiles)
            fps = AllChem.RDKFingerprint(mol, fpSize=bit_size)
            fpt = [int(x) for x in list(fps.ToBitString())]                               
        else:
            fpt = numpy.zeros(1025, dtype=numpy.int)
            # sm = [pybel.readstring("smi", x).write("can") for x in [smiles]]
            smcan = [pybel.readstring("smi", x) for x in [smiles]]
            fps = [x.calcfp(fptype=fptt) for x in smcan]
            fpt[[fps[0].bits]] = 1
            fpt = numpy.delete(fpt, 0)
            fpt = fpt.tolist()
        return fpt
    except:
        print("Oops!", sys.exc_info()[0], "occured.")
        print("Next entry.")
        print()
        #fpt = numpy.full((bit_size,1), 0)
        fpt = numpy.zeros(bit_size, dtype=numpy.int)
        return fpt

### function to convert SMILES to inchikeys
def smiles2inchikey(smiles):
    try:
        mol = AllChem.MolFromSmiles(smiles)
        inchikey = AllChem.MolToInchiKey(mol)
        return inchikey           
    except:
        print("Error converting SMILES to inchi", sys.exc_info()[0])
        return "-666"

def inchify_smiles(smiles):
    try:
        # Openbabel version, output should match obabel -:'<SMILES>' -osmi -xI
        mol = pybel.readstring('smi', smiles)
        inchi = mol.write('inchi').rstrip()
        inchikey = mol.write('inchikey').rstrip()
        inchified_smiles = pybel.readstring('inchi', inchi).write('smi').rstrip()
        
        # RDKIT equivalent, note the smiles output might differ from obabel
        #mol = AllChem.MolFromSmiles(smiles)
        #inchi = AllChem.MolToInchi(mol)
        #inchikey = AllChem.MolToInchiKey(mol)
        #mol2 = AllChem.MolFromInchi(inchi)                  
        #inchified_smiles = AllChem.MolToSmiles(mol2)        
        return {'inchified_smiles' : inchified_smiles,
                'inchi' : inchi,
                'inchikey' : inchikey}
    except:
        print("Error converting SMILES to inchi", sys.exc_info()[0])
        return {'inchified_smiles' : '-666',
                'inchi' : '-666',
                'inchikey' : '-666'}
    
### function to convert SMILES to binary fingerprints (uses RDKit)
def smiles2fpt_fcfp(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    fps = AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=1024,useFeatures=True)
    fpt = [int(x) for x in list(fps.ToBitString())]
    return fpt

### function to convert SMILES to binary fingerprints (uses RDKit)
def smiles2fpt_ecfp(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    fps = AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=1024,useFeatures=False)
    fpt = [int(x) for x in list(fps.ToBitString())]
    return fpt

### function to read in a CSV file
def read_csv_generic(fname, delimiter):
    idlist = []
    with io.open(fname, 'r', encoding='utf-8') as f:
        smil = csv.reader(f, delimiter=delimiter)
        for row in smil:
            idlist.append(row[0])
    return idlist



def unicode_csv_reader(unicode_csv_data, dialect=csv.excel, **kwargs):
    # csv.py doesn't do Unicode; encode temporarily as UTF-8:
    csv_reader = csv.reader(utf_8_encoder(unicode_csv_data),
                            dialect=dialect, **kwargs)
    for row in csv_reader:
        # decode UTF-8 back to Unicode, cell by cell:
        yield [unicode(cell, 'utf-8') for cell in row]

def utf_8_encoder(unicode_csv_data):
    for line in unicode_csv_data:
        yield line.encode('utf-8')

def csv_reader(csv_data, dialect=csv.excel, **kwargs):
    csv_reader = csv.reader(csv_data,
                            dialect=dialect, **kwargs)
    for row in csv_reader:
        yield row

### function to read in a CSV file
def read_csv_smiles(fname, delimiter):
    idlist = []
    with io.open(fname, 'r') as f:
        smil = csv_reader(f, delimiter=delimiter)
        keys = next(smil)
        if 'canonical_smiles' not in keys:
            print('There is no field with \'canonical_smiles\' in the header of your CSV file')
            sys.exit()
        nkeys = len(keys)
        rcount = 0
        for row in smil:
            rcount += 1
            if nkeys != len(row):
                print('Number of elements in row {0:d} does not much the number of keys'.format(rcount))
                print('Skipping this row')
            else:
                d = {k: v for e in keys for k, v in zip(keys, row)}
            idlist.append(d)
    return idlist


### Reads in a tab delimited file with pert_id and a fingerprint in hex format
def read_csv(fname, ftype, delimiter):
    olist = []
    indl = []
    cl = 0
    with open(fname, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            cl += 1
            if cl == 1:
                nf = len(row)
                indl = row
                ind = (i for i, val in enumerate(indl) if match('finger', val))
                for i in ind:
                    indl[i] = 'binary_fpt_old'
                print(indl)
            else:
                odict = {}
                for ii in range(nf):
                    if indl[ii] == 'binary_fpt_old':
                        odict[str(indl[ii])] = [int(x) for x in list(row[ii])]
                    else:
                        odict[str(indl[ii])] = row[ii]
                olist.append(odict)
    return olist


### Writes a tab delimited file with pert_id and binary fpt
def write_csv(binfpt, fname, delimiter):
    with open(fname, 'w') as f:
        binfptwriter = csv.writer(f, delimiter=delimiter)
        binfptwriter.writerow(binfpt)


### Writes a tab delimited file with pert_id and binary fpt
def write_header2csv(mydict, fname, delimiter):
    with open(fname, 'w') as f:
        binfptwriter = csv.writer(f, delimiter=delimiter)
        binfptwriter.writerow(mydict)


### Writes a tab delimited file with pert_id and binary fpt
def write_dict2csv(mydict, fname, delimiter):
    with open(fname, 'a') as f:
        binfptwriter = csv.writer(f, delimiter=delimiter)
        binfptwriter.writerow(mydict)


### Converts hex to bin
def hex2bin(hexstr):
    binstr = []
    nobits = 4
    hexstr = list(hexstr)

    for i in range(len(hexstr)):
        binstr.append(bin(int(hexstr[i], 16))[2:].zfill(nobits))
    binstr = ''.join(binstr)
    return binstr


### Read in a file with
def csv2binfpt(fname):
    hexfpt = read_csv(fname, '\t')
    binfpt = [[x, hex2bin(hexfpt[x])] for x in iter(hexfpt)]
    outname = fname[:-4] + '_binary.fpt'
    write_csv(binfpt, outname, '\n')
    return binfpt


### function to compare two lists (list1 is a list of dictionaries, like in dataJSON)
def compare_lists(list1, list2, key):
    reduce_list1 = []
    for i in range(len(list1)):
        if list1[i][key] in list2:
            print("{0:s} yes".format(i))
            reduce_list1.append(list1[i])
        else:
            print("{0:s} no".format(i))
    return reduce_list1


### Function to construct GCT structure and write gctx file in python
def write_gctx(data, ofile):
    dataJON = {}
    # Create a numpy matrix
    bfpt = [data[x]['binary_fpt'] for x in range(len(data))]
    bfpta = numpy.transpose(numpy.array(bfpt, dtype = 'i8'))
    print(bfpta.shape)
    #bfpta = numpy.transpose(numpy.array(bfpt, dtype='bool'))
    # Create column desc
    dataJO = copy.deepcopy(data)
    #[dataJO[x].pop('pert_id') for x in range(len(dataJO))]
    #[dataJO[x].pop('binary_fpt') for x in range(len(dataJO))]
    #for dkey in dataJO[0].keys():
    dataJON['pert_iname'] = [dataJO[x]['pert_iname'] for x in range(len(dataJO))]
    # Create cid and rid
    cid = [data[x]['pert_id'] for x in range(bfpta.shape[1])]
    rid = ["bit" + str(x + 1) for x in range(bfpta.shape[0])]
    data_df = pd.DataFrame(bfpta, 
                index = pd.Index(rid, name = "rid"),
                columns = pd.Index(cid, name = "cid"))
    #print(data_df.head())
    # TOADD Column metadata df
    row_df = pd.DataFrame(index = rid)
    col_df = pd.DataFrame(index = cid)

    #print(row_df)
    #print(col_df)

    #gcto(bfpta, rid, cid, {}, dataJON)
    #gcto.write(ofile, 'gctx')
    gco = gctoo.GCToo(data_df = data_df, row_metadata_df=row_df, col_metadata_df=col_df)
    #wg.write(gco, ofile)
    wgx.write(gco, ofile)

### Function to write output file in json, csv, or gctx format
def write_output(data=None, pert_id=None, ifile=None, ojson=None, ocsv=None, ogctx=None):
    if ojson is not None:
        write_json(data, ojson)
        print('Output written to JSON file: {0:s}'.format(ojson))
    elif ocsv is not None:
        write_header2csv(data[0].keys(), ocsv, '\t')
        for i in range(len(data)):
            write_dict2csv(data[i].values(), ocsv, '\t')
        print('Output written to CSV file: {0:s}'.format(ocsv))
    elif ogctx is not None:
        write_gctx(data, ogctx)
        print('Output written to GCTX file: {0:s}'.format(ogctx))
    else:
        if ifile is None:
            ogctx = pert_id + '_fpt.gctx'
            write_gctx(data, ogctx)
            print('Output written to GCTX file: {0:s}'.format(ogctx))
        else:
            ogctx = ifile.replace('.csv', '').replace('.json', '') + '_fpt.gctx'
            write_gctx(data, ogctx)
            print('Output written to GCTX file: {0:s}'.format(ogctx))
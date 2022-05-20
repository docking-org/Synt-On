from cmath import log
from xmlrpc.server import SimpleXMLRPCServer
import datetime, re
from rdkit import Chem
from rdkit.Chem import AddHs, AllChem
from rdkit.Chem.rdmolops import *
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from SyntOn_BulkFragmentationEnumerationAndAnaloguesDesign import main
import datetime, os, re, sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
srcPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.insert(1, srcPath)
from src.UsefulFunctions import *
from src.SyntOn import *


class SyntOnServer:
    def __init__(self, lib):
        self.lib = lib
        self.SynthLibrary = self.__buildLib(lib)
        
    def getlib(self):
        return self.lib

    def analog(self, input, output):
        return main(inp=input, inp2=None, strictAvailabilityMode=False, SynthLibrary = self.SynthLibrary, outDir=output, simTh=0.5, nCores=10,
            analoguesLibGen= True, maxNumberOfReactionCentersPerFragment=1, enumerationMode=True,
            parsedSynthLib=True)

    def __buildLib(self,lib, Ro2Filtration = False, findAnaloguesOfMissingSynthons = True):
        print("Processing BB library. It may take a few minutes, depending on the library size")
        fragBegTime = datetime.datetime.now()
        availableSynthons = {}
        pat = re.compile("\[\w*:\w*\]")
        for line in open(lib):
            sline = line.strip()
            if sline:
                mol = Chem.MolFromSmiles(sline.split()[0])
                if Ro2Filtration:
                    mol = AddHs(mol)
                    MolW = ExactMolWt(mol)
                    LogP = MolLogP(mol)
                    HDC = CalcNumHBD(mol)
                    HAC = CalcNumHBA(mol)
                    if MolW>200 or LogP > 2 or HDC > 2 or HAC > 4:
                        continue
                availableSynthons[sline.split()[0]] = {}
                availableSynthons[sline.split()[0]]["BBs"] = sline.split()[1]
                if findAnaloguesOfMissingSynthons:
                    availableSynthons[sline.split()[0]]["fp_b"] = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                availableSynthons[sline.split()[0]]["n_atoms"] = mol.GetNumAtoms()
                availableSynthons[sline.split()[0]]["n_rings"] = CalcNumRings(mol)
                availableSynthons[sline.split()[0]]["marks"] = sorted(
                    [sline.split()[0][m.start():m.start() + 2] + sline.split()[0][m.end() - 4:m.end()] for m in
                        re.finditer(pat, sline.split()[0])])
                availableSynthons[sline.split()[0]]["marksVallences"] = "+".join(sorted([atom.GetSymbol() + ":" + 
                            str(atom.GetTotalDegree()) for atom in mol.GetAtoms() if atom.GetAtomMapNum() != 0]))
        print("Lib BB reading time:")
        print(datetime.datetime.now() - fragBegTime)
        return availableSynthons

if __name__ == "__main__":     
    server = SimpleXMLRPCServer(("localhost", int(sys.argv[1])), logRequests= True, allow_none=True)
    try:
        server.register_instance(SyntOnServer(lib=sys.argv[2]))
        print("server is listening on port " + sys.argv[1])
        server.serve_forever()
    except Exception as e:
        print(e)
        print("server error")
        
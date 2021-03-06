from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.HeppyCore.statistics.tree import Tree as Tree
from ROOT import TFile

class TreeAnalyzerNumpy( Analyzer ):
    """Base TreeAnalyzerNumpy, to create flat TTrees.

    Check out TestTreeAnalyzer for a concrete example.
    IMPORTANT: FOR NOW, CANNOT RUN SEVERAL TreeAnalyzers AT THE SAME TIME!
    Anyway, you want only one TTree, don't you?"""

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(TreeAnalyzerNumpy,self).__init__(cfg_ana, cfg_comp, looperName)



    def beginLoop(self, setup) :
        super(TreeAnalyzerNumpy, self).beginLoop(setup)
        print setup.services
        if "outputfile" in setup.services:
            print "Using outputfile given in setup.outputfile"
            self.file = setup.services["outputfile"].file
        else :
            fileName = '/'.join([self.dirName,
                             'tree.root'])
            isCompressed = self.cfg_ana.isCompressed if hasattr(self.cfg_ana,'isCompressed') else 1
            print 'Compression', isCompressed
            self.file = TFile( fileName, 'recreate', '', isCompressed )
        self.tree = Tree('tree', self.name)
        self.declareVariables(setup)
        
    def declareVariables(self,setup):
        print 'TreeAnalyzerNumpy.declareVariables : overload this function.'
        pass

    def write(self, setup):
        super(TreeAnalyzerNumpy, self).write(setup)
        if "outputfile" not in setup.services:
            self.file.Write() 


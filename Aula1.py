from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta.teaching import *
import time as time_module
from pyrosetta import PyMOLMover
import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from rosetta.protocols import minimization_packing as pack_min
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.protocols.scoring import Interface
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.relax import FastRelax
import pyrosetta.rosetta.protocols.rigid as rigid_moves
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.select import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from pyrosetta.bindings.energies import *
from functions import PDB_pose_dictionairy
import pandas as pd
import csv

pyrosetta.init()
######
###### Carregando PDB e inserindo em objeto POSE
pose_norelax = pose_from_pdb("7t9k_alfa1_rbd.pdb.pdb")
pose = pose_from_pdb("7t9k_alfa1_rbd_relaxed.pdb")

###### Utilizando CLEAN pose para tratar possíveis erros na estrutura ao
###### carregar a pose.
###cleanATOM("7t9k_alfa1_rbd_relaxed")
pose_clean = pose_from_pdb("7t9k_alfa1_rbd_relaxed.clean.pdb")


###### Atribuindo Sequência da pose a uma variável
seq_original = pose.sequence()
seq_clean = pose_clean.sequence()

seq_original
###### Comparando sequência clean x original
print("Sequência sofreu alterações?",seq_original != seq_clean)
#########################################################################################
###### Trabalhando com Pose Residues
pose.total_residue()


###### Pegando cadeia do residuo e index
pose.pdb_info().chain(1)

###### Pegando numeração do resíduo (pose -> PDB)
pose.pdb_info().number(1)

###### Pegando numeração do resíduo (PDB -> pose)

pose.pdb_info().pdb2pose('A', 331)
pose.pdb_info().pose2pdb(1)

###### Dataframe mostrando index pdb x pose
## Function to create a dictionary containing |Chain|PDB Residue index| Pose Residue Index
def PDB_pose_dictionairy(pose_to):
    vectorIndexChain = []
    vectorIndexPDB = []
    vectorIndexPose = []
    ## Creating dictionary for residue position - PDB and Pose numbering
    for i in range(pose_to.total_residue()):
        vectorIndexChain.append(pose_to.pdb_info().chain(i+1))
        vectorIndexPDB.append(pose_to.pdb_info().number(i+1))
        vectorIndexPose.append(pose_to.pdb_info().pdb2pose(vectorIndexChain[i], vectorIndexPDB[i]))

    ## Inserting values into a pandas dataframe

    df_pdbtopose_dictionary = {"Chain" : vectorIndexChain,
                              "IndexPDB" : vectorIndexPDB,
                              "IndexPose": vectorIndexPose}
    df_dictionary = pd.DataFrame(df_pdbtopose_dictionary)
    return(df_dictionary)

df = PDB_pose_dictionairy(pose)
df

##### Energias no PyRosetta
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

##### Printando energia da pose - Com x Sem fastrelax
scorefxn.score(pose)
scorefxn.score(pose_norelax)

##### Decompondo as contribuições energéticas em uma tabela
scorefxn.show(pose)
pose.energies().residue_total_energy(1)

###### Retirando contribuição de termos especificos para cada residuo
###### Obs.: devolve os scores sem o weight
pose.energies().residue_total_energies(1)[fa_atr]
pose.energies().residue_total_energies(1)[fa_rep]
pose.energies().residue_total_energies(1)[fa_sol]



###### Calculando energias entre pares de resíduos

##### Printando pares de interações dado um threshold
pyrosetta.toolbox.atom_pair_energy.print_residue_pair_energies(
    1, ### posição para avaliar pares de interação
 	pose, ### pose
 	scorefxn, ### score
 	fa_atr, ### Termo de energia para printar
 	1 ### Threshold | valor minimo de energia para printar
)		



##### Inserindo valores em uma variavel
sss = pyrosetta.toolbox.atom_pair_energy._reisude_pair_energies(
    1,
 	pose,
 	scorefxn,
 	fa_rep,
 	threshold = 0
)		

for i, k in sss:
    print(i, k)


#######
for i in range(1,5):
    sss = pyrosetta.toolbox.atom_pair_energy._reisude_pair_energies(
    i,
 	pose,
 	scorefxn,
 	fa_rep,
 	threshold = 0)
    for j,k in sss:
        print("Pair:", i,j, "\nfa_rep: ", k,"\n")





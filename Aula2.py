from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.protocols.simple_moves import SmallMover
from pyrosetta.rosetta.protocols.simple_moves import ShearMover
from rosetta.protocols import minimization_packing as pack_min
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
#Core Includes
from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *

#Protocol Includes
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax


#### Inicializando pyrosetta
pyrosetta.init()
###### Carregando PDB e inserindo em objeto POSE
pose = pose_from_pdb("./Inputs/7t9k_alfa1_rbd_relaxed.pdb")
clone_pose = Pose()
clone_pose.assign(pose)
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

#### Movers

#### Criando o Movemap
movemap = MoveMap()


#### Atribruindo diferentes parâmetros ao MoveMap
movemap.set_bb(2, True)
movemap.set_chi(2,True)
movemap.set_bb_true_range(5,10)
movemap.set_jump(3, True)
movemap.show(10)
#### SmallMover e ShearMover

#### SmallMover
#### SmallMover gera perturbações no φi and ψi,

movemap = MoveMap()
movemap.set_bb_true_range(201,239)
n_moves = 5
kT = 1.0
smallmover = SmallMover(movemap, kT, n_moves)
print(smallmover)
smallmover.angle_max("E", 5) #beta-strand residues
smallmover.angle_max("H", 30) #helix residues
smallmover.angle_max("L", 20) #loop residues

#### Printando parametros do smallmover
print(smallmover)
##### Aplicando smallmover a uma pose
smallmover.apply(clone_pose)
##### Comparando energias pose antes e depois de aplicar mover
scorefxn.score(pose)
scorefxn.score(clone_pose)
#clone_pose.dump_pdb("./testemover.pdb")

#### ShearMover
#### ShearMover gera perturbações no φi e ψi-1.
shearmover = ShearMover(movemap, kT, n_moves)
shearmover.angle_max(10)
shearmover.angle_max("E", 5) #beta-strand residues
shearmover.angle_max("H", 30) #helix residues
shearmover.angle_max("L", 20) #loop residues
#### Printando parametros do shearmover
print(shearmover)
#### Aplicando shearmover a pose
clone_pose = Pose()
clone_pose.assign(pose)
shearmover.apply(clone_pose)

#### Comparando energia antes e depois de aplicar shearmover

scorefxn.score(pose)
scorefxn.score(clone_pose)
#clone_pose.dump_pdb("./testeshearmover.pdb")

#### Repacker

clone_pose = Pose()
clone_pose.assign(pose)

#### TaskFactory usado para manipular o repacker

tf = TaskFactory()
tf.push_back(operation.InitializeFromCommandline())
#### Converte o taskfactory em packer para poder analisar
packer_task = tf.create_task_and_apply_taskoperations(clone_pose)
print(packer_task)

#### Função que restringe a função repacker apenas para packing, desabilita design

tf.push_back(operation.RestrictToRepacking())
packer_task = tf.create_task_and_apply_taskoperations(clone_pose)
print(packer_task)

#### Criando Mover do PackRotamers

packer = pack_min.PackRotamersMover()
packer.task_factory(tf)

#### Aplicando packer a pose
pose_repacked = pose_from_pdb("./Inputs/6lzg.pdb")
before_energy = scorefxn.score(pose_repacked)
packer.apply(pose_repacked)

#### Comparando energias antes e após repacking
after_energy = scorefxn.score(pose_repacked)

print("Energy Value before repacking: {} REU\nEnergy Value after repacking:  {} REU".format(before_energy, after_energy))

##############

#### Packing regional
#### Utilização de ReisdueIndexSelector, para seleciona residuos de interesse
selected_region = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
selected_region.set_index("1-5")


#### Selecionando os residuos vizinhos aos que foram selecionados para também realizar o repack
#### Obs.: distância default utilizada para vizinhos é de 10A
#### https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.core.select.residue_selector.html

nbr_selector = selections.NeighborhoodResidueSelector()
nbr_selector.set_focus_selector(selected_region)
nbr_selector.set_include_focus_in_subset(True)

prevent_repacking_rlt = operation.PreventRepackingRLT()
prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
#tf.push_back(operation.RestrictToRepacking())


#### Printando os residuos selecionados para fazer o repacking
selecionados = []
print("Selected")
for i in select.get_residue_set_from_subset(nbr_selector.apply(pose)):
    print(i)
    selecionados.append(i)



packer_task = tf.create_task_and_apply_taskoperations(pose)
print(packer_task)



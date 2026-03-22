from utils.print_utils import print_dict
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import import_ipynb


# In[ ]:





# In[2]:


from openfile import open_log
from openfile import open_xyz


# In[3]:




# In[4]:


lines=open_log("fb_pyrrole.log",encoding="utf-8")
rows,data=open_xyz("pyrrole.xyz")


# In[5]:


import collections
#file = open_log()
#lines=file.readlines()
#file=open(filename,"r"
## get the M
def MO_coeff():
    MO_dict = collections.defaultdict(dict)
    MO_list = []
    MO_level= 0
    start=0

    for line in lines:
 
        if ('MO') in line:
            #if start==1:
            #    print("\nmo level: ", MO_level)
            #    print(MO_dict[MO_level])
     
            if ('Occ: 2') in line:
                line = line.strip().split()
                MO_level=int(line[2])
                MO_list.append(MO_level)
                start=1
            else:
                start=0
            
            
        elif start==1:
            line = line.strip().split()
            coeff = line[-1]
            atom_label = ("-").join(line[0:2])
            #atom_label=("-").join(line[2:0])
            
            if MO_level not in MO_dict:
                MO_dict[MO_level][atom_label]= float(coeff)
            elif atom_label not in MO_dict[MO_level]:
                MO_dict[MO_level][atom_label]= float(coeff)
            else:
                MO_dict[MO_level][atom_label]= MO_dict[MO_level][atom_label] + float(coeff)
                
    return MO_dict, MO_list


# In[6]:


MO_dict, MO_list = MO_coeff()
print(MO_list)


# In[7]:


print(rows[1])
smi=rows[1].split()[1]
print(smi)


# In[8]:


import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Atom
#name=input('enter the smiles format: ')
def val_elec():
   # chcl3=Chem.MolFromSmiles(name)
    m=Chem.MolFromSmiles(smi)
   # m=Chem.MolFromSmiles('F[Si](F)(F)F')
    #m=Chem.MolFromSmiles('FS(F)(F)(F)(F)F') ## PCL5
    #m=Chem.MolFromSmiles('C(Cl)(Cl)(Cl)Cl')  ## C2F2H2-cis
    #m=Chem.MolFromSmiles('C=C=C')  ## allene
    #m=Chem.MolFromSmiles('C[N+](=O)[O-]') ## Nitromethane
    #m=Chem.MolFromSmiles('[O-][S++]([O-])(C)C') ## (CH3)2SO
    valm=Chem.Descriptors.NumValenceElectrons(m)
    valshell=valm//2
    #print(valchcl3//2)
    return valshell


# In[9]:


#val_elec()


# In[10]:


def modified_mo_list():
    global MO_list
    val=val_elec()
    #val=9
    print("Valency: ", val)
    MO_list = MO_list[-val:]
    
modified_mo_list()
print("Modified mo list ", MO_list)


# In[11]:


all_atom_dict = collections.defaultdict(list)
mo_index_all = collections.defaultdict(list)

def atom_types():
    number=0
    H={'H':0}
    first_row={'C':0,'O':0,'F':0,'N':0,'B':0}
    second_row={'Na':0,'Al':0,'Si':0,'P':0,'S':0,'Cl':0}
    
    # one MO will have all the atom types hence, MO_dict[0]
    for atom in MO_dict[0].keys():
        
        # Initialize all atoms for (MO, coeff) list
        if atom not in all_atom_dict:
            all_atom_dict[atom] = []
            mo_index_all[atom] = []
          
        # Find no. of atoms
        atom = atom.split('-')[1]
        #print(atom)
        if atom in H:
            H[atom] += 1
        elif atom in first_row:
            first_row[atom] += 1
        else:
            second_row[atom] += 1
        
    return H, first_row, second_row        
    #print(H, first_row, second_row)


# In[12]:


atom_types()


# In[13]:


H=['H']
first_row=['C','N','O','F','B']
second_row=['Al','Si','P','S','Cl']


def create_mo_coeff_pair():
    global all_atom_dict
    
    for mo in MO_dict.keys():
        if mo not in MO_list:
            continue
        #print("\n", mo,)
        for atom in MO_dict[mo].keys():
            #print(atom)
            if atom in all_atom_dict.keys():
                all_atom_dict[atom].append((mo, MO_dict[mo][atom]))
                #mo_index_all[atom].append(mo)
                
                #print (atom, (mo, MO_dict[mo][atom]))


        
    


# In[14]:


create_mo_coeff_pair()


# In[15]:


#print_dict(all_atom_dict)


# In[16]:




def get_top_atom_coeff(central_atom, num_coeff,atom_type,x):
    top_atom_coeff_dict1 = collections.defaultdict(list)

    for atom in all_atom_dict:
        temp = atom.split('-')[1]
        
        if central_atom == atom:
            temp_list = get_top_n_ele(all_atom_dict[atom], num_coeff)
        elif atom_type == atom:
            temp_list = get_top_n_ele(all_atom_dict[atom], x)
        elif 'H' in atom:
            temp_list = get_top_n_ele(all_atom_dict[atom], 1)
        elif 'Al'==temp or 'B' == temp:
            temp_list = get_top_n_ele(all_atom_dict[atom], 3)
        else: ## for first and second row
            temp_list = get_top_n_ele(all_atom_dict[atom], 4)
        
        top_atom_coeff_dict1[atom] = temp_list
    return top_atom_coeff_dict1


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





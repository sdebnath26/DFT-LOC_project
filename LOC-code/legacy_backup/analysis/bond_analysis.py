from utils.print_utils import print_dict
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import import_ipynb
import collections
from collections import defaultdict
from collections import Counter
import itertools
import functools
import operator
import rdkit
from rdkit import Chem


# In[2]:


## function to create dictionary
dic_atom={}
dic_atomo_neigh={}
## Function to pretty prininting of dictinonary


# In[3]:


from openfile import open_xyz
from openfile import open_log
from Open_MO_all_atom_dict import get_top_atom_coeff
from Open_MO_all_atom_dict import all_atom_dict
from Open_MO_all_atom_dict import get_top_n_ele
from Open_MO_all_atom_dict import smi


# In[4]:


from geom_analysis import bond_distance 
from geom_analysis import angle_dict1
from geom_analysis import connectivity_dict
from geom_analysis import z_dict
from geom_analysis import sum


# In[5]:


#file_xyz=open_xyz("Propane-new.xyz")


# In[6]:


#print(file_xyz)


# In[7]:


top_atom_coeff_dict1 = get_top_atom_coeff(None,0,None,0)


# In[8]:


#print_dict(top_atom_coeff_dict1)


# In[9]:


def atom_index_pair():
    
    keys_index=[]
    keys_index_pair={}
    for keys in top_atom_coeff_dict1.keys():
    #print(top_atom_coeff_dict[keys])
        temp_index=[]
        for i in range(len(top_atom_coeff_dict1[keys])):
            #print(keys,top_atom_coeff_dict[keys][i][0])
            temp_index.append(top_atom_coeff_dict1[keys][i][0])
            #print(temp_index)
            keys_index_pair[keys]=temp_index
    return keys_index_pair


# In[10]:


print_dict(all_atom_dict)


# 

# In[11]:


print_dict(top_atom_coeff_dict1)


# In[12]:


print_dict(connectivity_dict)


# In[13]:


def index_match(top_atom_coeff_dict1):
    keys_index=[]
    keys_index_pair={}
    for keys in top_atom_coeff_dict1.keys():
    #print(top_atom_coeff_dict[keys])
        temp_index=[]
        for i in range(len(top_atom_coeff_dict1[keys])):
            #print(keys,top_atom_coeff_dict1[keys][i][0])
            temp_index.append(top_atom_coeff_dict1[keys][i][0])
            keys_index_pair[keys]=temp_index
    bd=[]
    bd =list(bond_distance.keys())   
    x=[]
    first_atom=[]
    second_atom=[]
    bond_pair_dict=collections.defaultdict(list)
    #print(bd)
    for i in range(len(bd)):
        first,second=bd[i].split(':')
        if first or second in keys_index_pair.keys():
            #print('First_atom:',first,keys_index_pair[first],'\t','Second_atom',second,keys_index_pair[second])
            first_atom.append(first)
            second_atom.append(second)
            x=keys_index_pair[first]
            #print(x)
            y=keys_index_pair[second]
            index_common=list(set(x).intersection(y))
            #print('bond_pair',bd[i],'\t','index_common:',index_common,'\n')
            bond_pair_dict[bd[i]]=index_common
            #print(bd[i],index_common)
            #return bd[i],index_common
    return bond_pair_dict,first_atom,second_atom


# In[14]:


index_match(top_atom_coeff_dict1)
bond_pair_dict={}
bond_pair_dict,first,second=index_match(top_atom_coeff_dict1)
#print(first)
print_dict(bond_pair_dict)


# In[15]:


## count the number of atoms 
bd=list(bond_distance.keys())

atom_count=defaultdict(list)
#print(first_atom,second_atom)
set1=first+second
#print(set1)
atom_count=collections.Counter(set1)

print_dict(atom_count)

central_atom = max(atom_count, key=atom_count.get)
print(central_atom)


# In[ ]:





# In[16]:


print_dict(bond_pair_dict)


# In[17]:


empty_key=[]
for keys in bond_pair_dict.keys():
    #print(keys,len(bond_pair_dict[keys]))#
    first1,second2=keys.split(':')
    if len(bond_pair_dict[keys])==0:
        #print(central_atom,keys)
        empty_key.append(keys)
        
print(empty_key)


# In[18]:


n_pair=[number for number in empty_key if empty_key.count(number)]
#print(len(n_pair))

update_bond_pair={}
## update the top_atom_coeff_dict
t_list=[]
key={}
for i in range(len(empty_key)):
   
    if central_atom in empty_key[i].split(':')[0] or empty_key[i].split(':')[1]:
#        #print(central_atom,all_atom_dict[central_atom])
        if(len(n_pair))==1:
            t_list=get_top_n_ele(all_atom_dict[central_atom], 5)
            top_atom_coeff_dict1[central_atom]=t_list
            print_dict(top_atom_coeff_dict1)
            x,first,second=index_match(top_atom_coeff_dict1)
            print('X',x)
#            #x={key:value}
            bond_pair_dict.update(x)
#            #print_dict(top_atom_coeff_dict1)
        elif (len(n_pair))==2:
            t_list=get_top_n_ele(all_atom_dict[central_atom],6)
            top_atom_coeff_dict1=t_list
            x,first,second=index_match(top_atom_coeff_dict1)
##            #x={key:value}
            bond_pair_dict.update(x)
            #key,value=index_match(top_atom_coeff_dict1)
            #x={key:value}
            #bond_pair_dict.update(x)
        else:
            print("No need to update MO ")


# In[19]:


bd_pair=[]
def bond_pair():
    for keys in bond_pair_dict.keys():
        keys
        x=len(bond_pair_dict[keys])
        a,b=keys.split(':')
        bd_p=list(zip((a,),(b,)))
        bd_pair.extend(bd_p*x)
    return bd_pair


# In[20]:


bond_pair()


# In[21]:


lp={}
list1=list(bond_pair_dict.keys())
print(list1)
c=[]
m=[]


# In[22]:


def bond_index1(bond_pair_dict):
    b_index1=defaultdict(list)
    first_atom=[]
    second_atom=[]
    x=[]
    y=[]
    c=[]
    temp_a=defaultdict(list)
    temp_b=defaultdict(list)
    for keys in bond_pair_dict.keys():
        a,b=keys.split(':')
        first_atom.append(a)
        second_atom.append(b)
        #print(a,b)
    #print(first_atom,second_atom)
        c.append(bond_pair_dict[keys])
    merged = list(itertools.chain.from_iterable(c))
    #print("merged",merged)
    pair1=list(zip(first,c))
    pair2=list(zip(second,c))
    final_pair=pair1+pair2
    #print("final",final_pair)
    for atom,index in final_pair:
        #print(atom,index)
        try:
            b_index1[atom].append(index)
        except KeyError:
            b_index1[atom]=index
    #print_dict(b_index1)
    flat_list=[]
    for keys in b_index1.keys():
        #print(keys,b_index1[keys])
        flat_list=[item for elem in b_index1[keys] for item in elem]
        #print("mod",keys,flat_list)
        x={keys:flat_list}
        #print(keys)
        b_index1.update(x)
    return b_index1


# In[23]:


bond_index1(bond_pair_dict)


# In[ ]:





# In[24]:


## test for duplicate pairs

top_atom_coeff_dict2_test={}
top_atom_coeff_dict2_test=top_atom_coeff_dict1.copy()


rev_dict={}
duplicates={}
dup_lis_1={}

list_a=[]
t_list=[]


for key, value in bond_pair_dict.items(): 
        for i in value:
            rev_dict.setdefault(i, set()).add(key)

result = [key for key, values in rev_dict.items() if len(values) > 1]


## print the duplicates

for key,values in rev_dict.items():
        if len(values) > 1:
            duplicates[key]=values
print_dict(duplicates)


for key, value in duplicates.items():
    list_a=[]
    list_dup_mod=[]
    for i in value:
        a,b=i.split(':')[0],i.split(':')[1]
        list_a.append(a)
        list_a.append(b)
        dup_lis_1[key]=list_a
        
    c=Counter(list_a)
    c_mod=c.copy()
    #print(c_mod)
    
## list_containing duplicates other than the central atom
list_dup_mod=[]
for key1 in c.keys():
        if c[key1] ==2:
            #print(key,c[key])
            c_mod.pop(key1)
            list_dup_mod=list(c_mod.keys())
            #print(list_dup_mod)
print(list_dup_mod) 
new_list=[]
for x in list_dup_mod:
    if "H" in x:
        new_list.append(x)
        
print(new_list)

## print the list containing H and non-hydrogen

#for i in range(len(list_dup_mod)):
#    print(list_dup_mod[i].split('-')[1])




# ###### py_list = ['a-1','b-2','c-3','a-4']
# new_list =[]
# for x in py_list: 
#     if "a" in x:
#        new_list.append(x)
# print(new_list)
# 

# In[25]:


top_atom_coeff_dict2={}
top_atom_coeff_dict2=top_atom_coeff_dict1.copy()


def pair_duplicates():
    
    rev_dict={}
    duplicates={}
    dup_lis_1={}


    list_a=[]
    t_list=[]
    for key, value in bond_pair_dict.items(): 
            for i in value:
                rev_dict.setdefault(i, set()).add(key)

    result = [key for key, values in rev_dict.items() if len(values) > 1]

    print(result)


    for key,values in rev_dict.items():
            if len(values) > 1:
                duplicates[key]=values
    print_dict(duplicates)

    for key, value in duplicates.items():
        list_a=[]
        list_dup_mod=[]
        for i in value:
            a,b=i.split(':')[0],i.split(':')[1]
            list_a.append(a)
            list_a.append(b)
            dup_lis_1[key]=list_a
        print(dup_lis_1)
        
        c=Counter(list_a)
        c_mod=c.copy()
        #print("Counter",c_mod)
        
        
        for key1 in c.keys():
            if c[key1] ==2:
            #print(key,c[key])
                c_mod.pop(key1)
                list_dup_mod=list(c_mod.keys())
        
        print(list_dup_mod)
        
        for x in range(len(list_dup_mod)):
            atom=list_dup_mod[x]                
            
            if list_dup_mod[0].split('-')[1]=="H" or list_dup_mod[1].split('-')[1]=="H":
                x=get_top_atom_coeff(None,0,atom,2)
                y={atom:x[atom]}
                top_atom_coeff_dict1.update(y)
            
            else:
                x=get_top_atom_coeff(None,0,atom,5)
                y={atom:x[atom]}
                top_atom_coeff_dict1.update(y)
       
        print_dict(top_atom_coeff_dict1)

    for i in range(len(t_list)):
        tup=[]
        print("t_list",t_list[i])
        for j in range(len(top_atom_coeff_dict1[t_list[i]])):
            if (t_list[i].split('-')[1]=="H") and top_atom_coeff_dict1[t_list[i]][j][0] not in result:
                w=top_atom_coeff_dict1[t_list[i]][j]
                tup.append(w)
                top_atom_coeff_dict2[t_list[i]]=tup
        
            elif (t_list[i].split('-')[1]=="C") and top_atom_coeff_dict1[t_list[i]][j][0] not in result:
                w=top_atom_coeff_dict1[t_list[i]][j]
                tup.append(w)
            #print("tup",tup)
                top_atom_coeff_dict2[t_list[i]]=tup
        
            z={t_list[i]:tup}
        print("z-dict")
        print_dict(z)
            
    return len(duplicates), top_atom_coeff_dict2


# In[26]:


pair_duplicates()


# In[27]:


#print(duplicates.values())


# In[ ]:





# In[28]:


## delete the item which appreared twice in the list
## pick one H first and 


# In[ ]:





# In[ ]:





# In[29]:


print_dict(top_atom_coeff_dict2)


# In[30]:


## create two duplicate num list a) one for H and other for non hydrogen
top_atom_coeff_dict2={}
top_atom_coeff_dict2=top_atom_coeff_dict1.copy()
def pair_duplicates_test():
    
    rev_dict={}
    duplicates={}

    for key, value in bond_pair_dict.items(): 
            for i in value:
                rev_dict.setdefault(i, set()).add(key)

    result = [key for key, values in rev_dict.items() if len(values) > 1]



    for key,values in rev_dict.items():
            if len(values) > 1:
                duplicates[key]=values
    print_dict(duplicates)

    if len(duplicates)!=0:
        temp_H=[]
        temp_other=[]
        dup_list_H=[]
        dup_list_other=[]
        dup_list_all=[]
    
        for keys,value in duplicates.items():
        
            #print(value)
                value1=list(value)
                temp_list=[]
                for i in range(len(value1)):
                    m=value1[i].split(":")
                    temp_list2=[]
                    c=keys
                    a,b=value1[i].split(':')
        #print(x)
                    temp_list2.append(a)
                    temp_list2.append(b)
                    temp_list.append(temp_list2)
           
            
                    for j in range(len(m)):
                        if "H" in m[j].split("-")[1]:
                        #print(keys,duplicates[keys])
                            temp_H.append(keys)
                            dup_list_H=list(set(temp_H))
                dup_list_all.append(temp_list)
    
        dup_list_other=list(set(dup_list_H).symmetric_difference(set(result)))
        print(result)
    
        t_list=[]
        
        for i in range(len(dup_list_all)):
            print(dup_list_all[i])
        
            if len(dup_list_all[i])==2:
                        check = any(item in dup_list_all[i][0] for item in dup_list_all[i][1])
                        t_list.append(dup_list_all[i][1][1])
            for i in range(len(t_list)):
                tup=[]

                if t_list[i].split('-')[1]=='H':
                    atom1=t_list[i]
                    x=get_top_atom_coeff(None,0,atom1,2)
        
                    y={atom1:x[atom1]}
            
                    top_atom_coeff_dict1.update(y)
                else:
                    atom1=t_list[i]
                    x=get_top_atom_coeff(None,0,atom1,5)
       
                    y={atom1:x[atom1]}
        
                    top_atom_coeff_dict1.update(y)
                    print_dict(top_atom_coeff_dict1)

            for j in range(len(top_atom_coeff_dict1[t_list[i]])):
                    print(t_list[i],top_atom_coeff_dict1[t_list[i]][j][0])
                    if (t_list[i].split('-')[1]=="H") and (top_atom_coeff_dict1[t_list[i]][j][0] not in dup_list_H):
                    #print("t",t_list[i],top_atom_coeff_dict1[t_list[i]][j][0])
                        w=top_atom_coeff_dict1[t_list[i]][j]
                        tup.append(w)
                        top_atom_coeff_dict2[t_list[i]]=tup
                    
                    elif  (top_atom_coeff_dict1[t_list[i]][j][0] not in result):
                    #print("other",t_list[i],top_atom_coeff_dict1[t_list[i]][j][0])
                        w=top_atom_coeff_dict1[t_list[i]][j]
                        tup.append(w)
                        print("tup",tup)
                        top_atom_coeff_dict2[t_list[i]]=tup
        
            z={t_list[i]:tup}
            print_dict(z)
        return len(duplicates), top_atom_coeff_dict2
    
    else:
        return len(duplicates)
            
        #top_atom_coeff_dict1.update(z)


# In[31]:


#top_atom_coeff_dict1.update(z)
#pair_duplicates()


# In[32]:


pair_duplicates_test()


# In[33]:


### Part from New index matching
#top_atom_coeff_test_2={}
#top_atom_coeff_dict2=top_atom_coeff_dict1.copy()
#def pair_duplicates():
    
#    rev_dict = {}
#    duplicates={}
# 
#    for key, value in bond_pair_dict.items(): 
#        for i in value:
#            rev_dict.setdefault(i, set()).add(key)
#    
#
#    result = [key for key, values in rev_dict.items() if len(values) > 1] #

#    for key,values in rev_dict.items():
#        if len(values) > 1:
#            duplicates[key]=values
#    print_dict(duplicates)
    
#    if len(duplicates)!=0:
#        
#        dup_list=[]
#        dup_num_list=[]
#        for keys,value in duplicates.items():
#            value1=list(value)
#            temp_list=[]
#   
#            for i in range(len(value1)):
#                temp_list2=[]
#                c=keys
#                a,b=value1[i].split(':')
#        #print(x)
#                temp_list2.append(a)
#                temp_list2.append(b)
#                temp_list.append(temp_list2)
#            dup_num_list.append(keys)
#            dup_list.append(temp_list)
#            print("dup_list",dup_list)
#            t_list=[]
#            for i in range(len(dup_list)):
#                print("dup_list",dup_list[i])
    #print(dup_list[i])
#                if len(dup_list[i])==2:
#                    check = any(item in dup_list[i][0] for item in dup_list[i][1])
#                    #print(dup_list[i][0][1],'\t',dup_list[i][1][1])
#                    t_list.append(dup_list[i][1][1])
#        for i in range(len(t_list)):
#            tup=[]
#    temp={}
#            if t_list[i].split('-')[1]=='H':
#                atom1=t_list[i]
#                x=get_top_atom_coeff(None,0,atom1,2)
#        #print(atom1,x[atom1])
#                y={atom1:x[atom1]}
#            #print(y)
#                top_atom_coeff_dict1.update(y)
#            else:
#                atom1=t_list[i]
#                x=get_top_atom_coeff(None,0,atom1,5)
        #print(atom1,x[atom1])
#                y={atom1:x[atom1]}
        #print(y)
#                top_atom_coeff_dict1.update(y)
#                print("dup_num_list",dup_num_list)
#            for j in range(len(top_atom_coeff_dict1[t_list[i]])):
#                if top_atom_coeff_dict1[t_list[i]][j][0] not in dup_list_H:
#                    print(t_list[i],":",top_atom_coeff_dict1[t_list[i]][j][0])
#                    w=top_atom_coeff_dict1[t_list[i]][j]
#                    tup.append(w)
#            z={t_list[i]:tup}
#            top_atom_coeff_dict1.update(z)
#        return len(duplicates), top_atom_coeff_dict1
#    else:
#        return len(duplicates)


# In[34]:


#pair_duplicates()


# In[ ]:





# In[35]:


print_dict(top_atom_coeff_dict2)


# In[36]:


print_dict(connectivity_dict)


# In[ ]:





# In[ ]:





# In[37]:


if pair_duplicates():
    print("correction needed")
    #print_dict(test_list)
    n,test_list=pair_duplicates()
    print_dict(test_list)
    bond_pair_dict1,first1,second1=index_match(test_list)
    print_dict(bond_pair_dict1)
    bd_pair.clear()
    bond_pair()
    b_index2=bond_index1(bond_pair_dict1)
    
else:
    bond_pair_dict1,first1,second1=index_match(top_atom_coeff_dict1)
    print_dict(bond_pair_dict1)
    bd_pair.clear()
    bond_pair()
    #b_index1.clear()
    b_index2=bond_index1(bond_pair_dict1)
    #continue
    #print_dict(top_atom_coeff_dict1)
    print("correction not needed")


# In[38]:


lp_dict=defaultdict(list)
def lone_pair():
    for atom in top_atom_coeff_dict1.keys():
        if pair_duplicates():
            n,test_list1=pair_duplicates()
            lp = [x[0] for x in test_list1[atom] if x[0] not in b_index2[atom]]
            lp_dict[atom]=lp
        else:
            lp = [x[0] for x in top_atom_coeff_dict1[atom] if x[0] not in b_index2[atom]]
            lp_dict[atom]=lp
    return lp_dict


# In[64]:


lone_pair()
#print_dict(b_index2)
print("lp")
print_dict(lp_dict)


# In[40]:


atom_list={'H':1,'C':4,'N':5,'O':6,'F':7,'Al':3,'Si':4,'P':5,'S':6,'Cl':7}
atom_net_charge={}

## Calculate Net_charge of each atom, if there is a bp then -1 contribution and if lp then -2 contribution
def net_charge_atom():
    for atom in top_atom_coeff_dict1.keys():
        #print(atom,len(bp_dict[atom]))
       # print("atom",atom,'\t',"lone-pair",len(lp_dict[atom]),'\t',"bond-pair",len(bp_dict[atom]))
        temp_label=atom.split('-')[1]
       # print(temp_label)
        n_lp=len(lp_dict[atom])
        n_bp=len(b_index2[atom])
        total_contribution=(-1)*n_bp+(-2)*n_lp
        if temp_label in atom_list:
            val=atom_list[temp_label]
            #print(val)
        n_charge=val+total_contribution
        print(atom,n_charge)
        atom_net_charge[atom]=n_charge
net_charge_atom()


# In[41]:


atom_valency={}
def valency_atom():
    for atom in top_atom_coeff_dict1.keys():
        temp_label=atom.split('-')[1]
        atom_val=abs(-2*((len(lp_dict[atom]))+len(b_index2[atom])))
        print(atom.split('-')[1],atom_val)
        atom_valency[atom]=atom_val
        if temp_label=='H' and atom_val==2:
            print('valency satisfied for H')
        elif temp_label !='H' and atom_val==8:
            print('octet rule filled')
        else:
            print('octet rule is not satisfied')
valency_atom()


# In[42]:


ct={}
ct_list=[]
def CT_atom_pair():
    for keys in bond_pair_dict.keys():
        a,b=keys.split(':')
        net=int(atom_net_charge[a])+int(atom_net_charge[b])
        print(a,b,atom_net_charge[a],atom_net_charge[b])
        
        #print(a,atom_net_charge[a],'\t',b,atom_net_charge[b])
        if int(atom_net_charge[a])==0 and int(atom_net_charge[b])==0:
            #print("Zero charges")
            continue
        elif net !=0:
            continue
            
        print(a,b,"Equal and opposite sign ")
        ct[keys]="CT-pair"
        i=(a,b)
        ct_list.append(i)
        return True
    return False


# In[43]:


CT_atom_pair()
print(ct_list)


# In[44]:


neighbour_charge={}
central_charge={}
ct_pair1=[]
def charge_sum():
    sum=0
    for atom in top_atom_coeff_dict1.keys():
        if atom!=central_atom:
            if atom_net_charge[atom]!=0:
                #print(central_atom,atom_net_charge[central_atom],atom,atom_net_charge[atom])
                sum=sum+int(atom_net_charge[atom])
                neighbour_charge[atom]=atom_net_charge[atom]
                #print(atom)
    net_sum=sum+int(atom_net_charge[central_atom])
    central_charge[central_atom]=atom_net_charge[central_atom]
    for atom in neighbour_charge.keys():
        print(atom,central_atom)
        ct_pair1.append((atom,central_atom))
    return net_sum,ct_pair1


# In[45]:


net_sum,CT_pair=charge_sum()
print(CT_pair)


# In[46]:


def b_order():
    s_list=[]
    d_list=[]
    t_list=[]
    if pair_duplicates():
        for b in bond_pair_dict1.keys():
            if len(bond_pair_dict1[b])==1:
            #print(b,bond_distance[b],'Single-bond','\n')
                s_list.append(b)
            elif len(bond_pair_dict1[b])==2:
            #print(b,bond_distance[b],'Double-bond','\n')
                d_list.append(b)
            elif len(bond_pair_dict1[b])==3:
            #print(b,bond_distance[b],'Triple-bond','\n')
                t_list.append(b)
    else:
         for b in bond_pair_dict.keys():
            if len(bond_pair_dict[b])==1:
            #print(b,bond_distance[b],'Single-bond','\n')
                s_list.append(b)
            elif len(bond_pair_dict[b])==2:
            #print(b,bond_distance[b],'Double-bond','\n')
                d_list.append(b)
            elif len(bond_pair_dict[b])==3:
            #print(b,bond_distance[b],'Triple-bond','\n')
                t_list.append(b)
    return s_list,d_list,t_list


# In[47]:


b_order()


# In[48]:


bd_pair.clear()
numbers=bond_pair()
print(numbers)
single = [number for number in numbers if numbers.count(number) == 1]
double = [number for number in numbers if numbers.count(number) == 2]
print(double)


# In[49]:


db_type={}
                
check_c_double_bond()
print_dict(db_type)


# In[50]:


def polar_F():
    d_1=list(set(double))
    count=0
    #print(numbers, "\n",d_1)
    if not check_c_double_bond():
        return 
    
    for ele in numbers:
        #print('ele','\t',ele)
        if ('F' in ele[0]) and ('C' == ele[1].split('-')[1]): 
            t1 = [ele[1] for x in d_1 if ele[1]==x[0] or ele[1]==x[1]]
            if(t1):
                count=count+1
            
        elif ('F' in ele[1]) and ('C' == ele[0].split('-')[1]) :
            t2 = [ele[0] for x in d_1 if ele[0]==x[0] or ele[0]==x[1]]
            if t2:
                count=count+1
    #print(count)     
    return count
polar_F()


# In[51]:


db_type={}
    
    for keys in bond_pair_dict.keys():
        a,b=keys.split(':')
        count=0
        if (len(bond_pair_dict[keys])==2):
            if 'C' in a.split('-')[1] and 'C' in b.split('-')[1]:
                count=count+1
                db_type['C=C']=count
            elif ('C' in a.split('-')[1] and 'O' in b.split('-')[1]) or ('O' in a.split('-')[1] and 'C' in b.split('-')[1]):
                count=count+1
                db_type['C=O']=count
    return True
                
check_c_double_bond()
print_dict(db_type)


# In[52]:


def bond_with_two_double_bonds():
    d_1=list(set(double))
    #print(d_1)
    flat_list = [item for sublist in d_1 for item in sublist]
    #print(flat_list)
    duplicates = [item for item, count in collections.Counter(flat_list).items() if count ==2]
    #print(len(duplicates))
    if len(duplicates)==1:
        #print("atom has two double bonds")
        return True
bond_with_two_double_bonds()


# In[53]:


def bond_order1(atom, numbers):
    #numbers = [('4-Cl', '0-P'), ('0-P', '3-Cl'), ('0-P', '2-Cl'), ('0-P', '1-O'),('0-P', '1-O'),('0-P', '1-O'),('0-P', '2-Cl')]
        
    single = [number for number in numbers if numbers.count(number) == 1 and (number[0]==atom or number[1]==atom)]
    double = [number for number in numbers if numbers.count(number) == 2 and (number[0]==atom or number[1]==atom)]
    triple = [number for number in numbers if numbers.count(number) == 3 and (number[0]==atom or number[1]==atom)]
    unknown = [number for number in numbers if numbers.count(number) > 3 and (number[0]==atom or number[1]==atom)]
    #print("\nAtom:", atom, "\n\tsingle:",single,"\n\tnumber of single bonds",len(single))
    #print("\n\n\tdouble:",double,"\n\tnumber of double bonds",len(list(set(double))))
    #print("\n\n\ttriple:",triple,"\n\tnumber of triple bonds",len(list(set(triple))))
    
    return single,len(single), double, len(list(set(double))), triple, len(list(set(triple)))


# In[54]:


m=Chem.MolFromSmiles(smi)

atom_n_ring=[]

for atom in m.GetAtoms():
    print(atom.GetIdx())
    
    
def isRingAromatic(mol,bondRing):
    for id in bondRing:
        if not mol.GetBondWithIdx(id).GetIsAromatic():
            return False
        return True
    
n_rings=Chem.GetSSSR(m)
print(n_rings)
id=()
ri=m.GetRingInfo() ## get the ring info
print(list(ri.AtomRings()))
idx_atom=ri.AtomRings() ## print the atom in rings
idx_ring=[]
count=0
#print(m.GetAtomWithIdx(1).IsInRing())

num_ele_x=[ele for tupl in idx_atom for ele in tupl]

idx=[]
def isRing(mol,ring):
    for idx in ring:
        if not mol.GetAtomWithIdx(idx).IsInRing():
            return False
        return True


# In[55]:


from geom_analysis import sum


# In[56]:


print(round(sum))


# In[57]:


atom_hybrid = {}
atom_correction = {}
def hybridisation():
    
    for atom in top_atom_coeff_dict1.keys():
        #print(atom)
        ele = atom.split('-')[1]
        if "H" in ele:
            continue
            
        s,sb,d,db,t,tb = bond_order1(atom, numbers)
        print(s)
        lp = len(lp_dict[atom])
        print ("\nAtom: ", atom, "\n\tsp,db,tb:", sb, db,tb, "\n\tLp:", lp)
        
        total = sb+db+tb+lp
        
        if ("N" in ele or "P" in ele):
            if total >= 3 and atom_net_charge[atom] == 1:
                print("\tHybridization: Quarternary")
                atom_hybrid[atom] = "Quarternary"   
                atom_correction[atom]=6.04
                print("\tCorrection: ",atom_correction[atom])
        if ("Cl" in ele or "P" in ele or "S" in ele) and atom_valency[atom] > 8:
            print("\tOctet expansion")
            atom_hybrid[atom] = "Octet Expansion"
            atom_correction[atom]= 4.84
            print("\tCorrection: ",atom_correction[atom])
        elif(total == 4):
            if ("N" in ele) and (round(sum)==360):
                print ("\tHybridization: sp2, conjugated")
                atom_hybrid[atom] = "sp2"
                atom_correction[atom]=4.46
                print("\tCorrection: ",atom_correction[atom])
            elif ("N" in ele or "P" in ele) and atom_net_charge[atom] != 1:
                print ("\tHybridization: sp3")
                atom_hybrid[atom] = "sp3"
                atom_correction[atom]=3.68
                print("\tCorrection: ",atom_correction[atom])
            elif ("O" in ele):
                print ("\tHybridization: sp3")
                atom_hybrid[atom] = "sp3"
                atom_correction[atom]=1.78
                print("\tCorrection: ",atom_correction[atom])
        elif (total == 3):
            print ("\tHybridization: sp2")
            atom_hybrid[atom] = "sp2"
            if ("N" in ele or "P" in ele) and atom_net_charge[atom] !=1:
                atom_correction[atom]=4.46
                print("\tCorrection: ",atom_correction[atom])
            elif ("O" in ele):
                atom_correction[atom]=0.95
                print("\tCorrection: ",atom_correction[atom])
        elif (total == 2):
            print ("\tHybridization: sp")
            atom_hybrid[atom] = "sp"
            if ("Be" in ele):
                atom_correction[atom]=7.31
                print("\tCorrection: ",atom_correction[atom])


# In[58]:


hybridisation()


# In[59]:


b_type_cyclic={}
def bond_corrections_cyclic():
    b_correction=False
    for number in numbers:
        bond_num=numbers.count(number)
        atom_temp_list= list(zip((number[0].split('-')[1],),(number[1].split('-')[1],)))
        if isRingAromatic(m,num_ele_x):
            
            if bond_num==1:
                if "H" in atom_temp_list[0]:
                    continue
                elif "F" in atom_temp_list[0]:
                    continue
                else:
                    b_type_cyclic[number]='single-bond in aromatic ring'
                    b_correction=True
            elif bond_num==2:
                b_type_cyclic[number]='double-bond in aromatic ring'
                b_correction=True   
        elif isRing(m,num_ele_x):
            if len(num_ele_x)==3:
                b_type[number]="Short, three membered ring"
                b_correction=True
                bond_cor[number]= -1.41
            else:
                b_type[number]="Medium, Ring C-C bond correction"
                b_correction=True
                bond_cor[number]= -1.96


# In[60]:


b_type={}
bond_cor={}
a_list={}
second_row=['Na','Mg','Al','Si','P','S','Cl','Ne']
def bond_corrections_final(number):
    ct_list_keys1={}
    b_correction=False
   # for number in numbers:
    bond_num=numbers.count(number)
    atom_list=list(zip((number[1],),(number[0],)))
    print(atom_list,CT_pair)
    atom_temp_list= list(zip((number[0].split('-')[1],),(number[1].split('-')[1],)))
    #a_list=dict.fromkeys(atom_list)
    #ct_list_keys1=dict.fromkeys(CT_pair)
    #print(ct_list_keys1.keys())
    if bond_num==1:
        if atom_list[0] in CT_pair:
            b_type[number]="single-bond with CT-character"
            b_correction=True
            bond_cor[number]=-4.53
        #if atom_temp_list[0]==ct_pair:
        #    b_type[number]="single-bond with CT-character"
        #elif isRing(m,num_ele_x):
        #    if len(num_ele_x)==3:
        #        b_type[number]="Short, three membered ring"
        #        b_correction=True
        #        bond_cor[number]= -1.41
        #    elif len(num_ele_x)==4:
        #        b_type[number]="Medium, Ring C-C bond correction"
        #        b_correction=True
        #        bond_cor[number]= -1.96
        #    elif not 
                
        elif isRingAromatic(m,num_ele_x):
            pass
        
        elif isRing(m,num_ele_x) and len(num_ele_x)==3:
                b_type[number]="Three membered ring"
                b_correction=True
                bond_cor[number]=-1.41
                            
        elif (('Li','Li') in atom_temp_list) or (('B','Cl') in atom_temp_list or ('Cl','B') in atom_temp_list) or (('C','N') in atom_temp_list or ('N','C') in atom_temp_list) or (('C','O') in atom_temp_list or ('O','C') in atom_temp_list):
            b_type[number]="Short, SSBC correction"
            b_correction=True
            bond_cor[number]= -1.41
            
        elif (('C','C') in atom_temp_list) or (('C','Cl') in atom_temp_list or ('Cl','C') in atom_temp_list) or (('N','N') in atom_temp_list) or (('O','O') in atom_temp_list) or (('N','O')in atom_temp_list or ('O','N') in atom_temp_list) or (('F','F') in atom_temp_list) or (('C','O') in atom_temp_list or ('O','C') in atom_temp_list) or (('Na','Na') in atom_temp_list):
            b_type[number]="Medium,MSBC correction"
            bond_cor[number]= -1.96
            b_correction=True
            
        elif (('Si','C') in atom_temp_list or ('C','Si') in atom_temp_list) or (('S','C') in atom_temp_list or ('C','S') in atom_temp_list) or (('S','O') in atom_temp_list or ('O','S') in atom_temp_list):
            b_type[number]="Long"
            bond_cor[number]=-2.63
            b_correction=True 
            
            
        if "H" in atom_temp_list[0]:
                #print("A-H bond")
            if "O" in atom_temp_list[0] or "Cl" in atom_temp_list[0] or "F" in atom_temp_list[0]:
                b_type[number]="polar"
                bond_cor[number]=-1.42
                b_correction=True
                
            elif "S" in atom_temp_list:
                b_type[number]="intermediate"
                bond_cor[number]=0
                b_correction=True
                
            else:
                b_type[number]="non-polar"
                bond_cor[number]=0.37
                b_correction=True
                
        if "F" in atom_temp_list[0]:
            if polar_F():
                b_type[number]="non-polar"
                bond_cor[number]=0.88
                b_correction=True
                
            elif "Li" in atom_temp_list[0] or "B" in atom_temp_list[0] or "C" in atom_temp_list[0] or "O" in atom_temp_list[0] or "Cl" in atom_temp_list[0]:
                b_type[number]="polar"
                bond_cor[number]=-0.95
                b_correction=True
                
            elif "N" in atom_temp_list[0]:
                b_type[number]="intermediate"
                bond_cor[number]=0
                b_correction=True
                
            elif "S" in atom_temp_list[0] or "P" in atom_temp_list[0] or "Si" in atom_temp_list[0]:
                b_type[number]="Si-F or S-F or P-F bond"
                bond_cor[number]=-4.45
                b_correction=True
                
    elif bond_num==2:
        if isRingAromatic(m,num_ele_x):
            pass
        
        elif  bond_with_two_double_bonds():
            b_type[number]="bond connected with another double-bond"
            b_correction = True
            
        else:
            b_type[number]="double-bond"
            b_correction = True
            bond_cor[number]=-1.03
            
    elif bond_num==3:
        if ("C" in atom_temp_list[0] and "C" in atom_temp_list[0]) or ("N" in atom_temp_list[0] and "N" in atom_temp_list[0]) or ("P" in atom_temp_list[0] and "P" in atom_temp_list[0]):
                b_type[number]="Non-polar triple bond"
                b_correction= True
                bond_cor[number]=-1.68
                
        elif ("C" in atom_temp_list[0] and "N" in atom_temp_list[0]) or ("C" in atom_temp_list[0] and "O" in atom_temp_list[0]) or ("C" in atom_temp_list[0] and "S" in atom_temp_list[0]) or ("Si" in atom_temp_list[0] and "O" in atom_temp_list[0])  :
                b_type[number]="Polar triple bond"
                b_correction=True
                bond_cor[number]=0.88
                
    return b_correction


# In[61]:


def pair_second_row():
    for number in numbers:
        bond_num=numbers.count(number)
        atom_temp_list= list(zip((number[0].split('-')[1],),(number[1].split('-')[1],)))
        atom_wo_split= list(zip((number[0],),(number[1],)))
        if set(atom_temp_list[0]).issubset(second_row):
            print(atom_temp_list,"Long bond, pair of second row elements")
            b_type[number]="Long bond, pair of second row elements"
            bond_cor[number]=-2.63
pair_second_row()


# In[62]:


for number in numbers:
    bond_corrections_final(number)
print_dict(b_type)
print_dict(bond_cor)
bond_corrections_cyclic()
print_dict(b_type_cyclic)


# In[ ]:





# In[ ]:





# In[ ]:





# In[63]:


## ESBC correction
#print(set(set(numbers)))
#print(ct_list)
esbc_list=[]
ct_list_keys={}
ct_list_keys=dict.fromkeys(ct_list)
for number in numbers:
    atom_t_list=list(zip((number[0].split('-')[1],),(number[1].split('-')[1],)))
    #print(atom_t_list)
    bond_num=numbers.count(number)
    a_list=list(zip((number[1],),(number[0],)))
    print(a_list)
    #print(a_list)
    #print(a_list)
    #print(bond_num)
    if bond_num==1:
        if "H" in atom_t_list[0]:
            #print(number,"Bond contains H: no ESBC correction")
            continue
        elif "F" in atom_t_list[0]:
           # print("Bond contains F: no ESBC correction")
            continue
        elif a_list[0] in CT_pair:
            print(number,"Bond is a CT-pair, No ESBC correction")
            continue
        #elif ct_list_keys.keys() & a_list.keys():
        #    print(number,"Bond is a CT-pair, No ESBC correction")
        #    continue
        elif (((len(num_ele_x)==3)) or (len(num_ele_x)==4)):
            continue 
        else :
            print(number)
            esbc_list.append(number)

## create key value pair having the connectivity info from the esbc list
esbc_dict_list=defaultdict(list)

for tup in esbc_list:
    print(tup)
    esbc_dict_list[tup[0]].append(tup[1])
    esbc_dict_list[tup[1]].append(tup[0])

print_dict(esbc_dict_list)

net_esbc=0
for keys in esbc_dict_list.keys():
    esbs_corr=0
    #print(keys,len(esbc_dict_list[keys]))
    n=len(esbc_dict_list[keys])
    esbc_corr=n*(n-1)
    print(keys,n,esbc_corr)
    net_esbc=net_esbc+esbc_corr

print('\n','net-ESBC:',net_esbc)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





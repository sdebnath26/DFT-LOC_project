from utils.print_utils import print_dict
#!/usr/bin/env python
# coding: utf-8

# In[1]:


## import functions 
import collections
import math
import numpy as np
import constants
import import_ipynb
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


#print(constants.Constants.m_dict) ## print dictionary from the Constants class of constant module)
dic_atom={}
dic_atomo_neigh={}
## Function to pretty prininting of dictinonary


# In[3]:


from openfile import open_xyz
from openfile import open_log


# In[4]:


rows,data=open_xyz("pyrrole.xyz")


# In[5]:


print(rows)


# In[6]:


## open the xyz
#with open("Propane-new.xyz", "r") as file:
#    read= file.read().strip()
#    print(read)
#    rows=read.split("\n")
#    data=[rows[0],rows[1],rows[2]]
#    data=data+[row for row in rows[2:] if row.strip !=""]


# In[ ]:





# In[7]:


##number of element
n=len(rows[3:])
## create a dictionary (x_dict), which contains name, number of atoms,atoms_list, xyz coordinate
x_dict={"n":n,
       "name":rows[1].strip(),
       "atoms":[],
       "xyz":[]}
if n!= len(rows)-3:
    print("invalid data format")
for i in range(n):
    sep_rows=rows[i+3].split()
    row =[x for x in sep_rows if x.strip() !=""]
    x_dict["atoms"].append(row[0].strip())
    coord = [float(j) for j in row[1:]]
    x_dict["xyz"].append(coord)
x_dict["xyz"]=np.transpose(x_dict["xyz"])
if len(x_dict["xyz"])!=3:
    print("invalid xyz format")


# In[8]:


print_dict(x_dict)


# In[9]:


print(x_dict.values())
for keys in x_dict.keys():
    print(len(keys))


# In[10]:


## indexing the atoms
index1=[]
#z_dict.clear()
for i in range(len(x_dict["atoms"])):
    #print(str(x_dict["atoms"][i])+"-"+str(i))
    y="".join(str(i)+"-"+str(x_dict["atoms"][i]))
    print(y)
    index1.append(y)
print(index1)


# In[11]:


z_dict={"atoms":[]}
z_dict["atoms"]=index1
#print_dict(z_dict)


# In[12]:


print(x_dict["atoms"])


# In[13]:


## Function to print the distance
def distance(a,b):
    xyz=np.transpose(x_dict["xyz"])
    a_array=xyz[a]
    b_array=xyz[b]
    #print(a_array,b_array)
    d = sum((a_array - b_array)**2)*0.5
    d1 = np.linalg.norm(a_array-b_array)
    return d,d1


# In[14]:


#distance(0,3)


# In[15]:


l_unit_dict = {"ang": 1.0,
                   "m": constants.Constants.ang,
                   "cm": constants.Constants.ang*100.0,
                   "a": constants.Constants.ang/constants.Constants.a}
m_unit_dict = {"u": 1.0,
                "kg": constants.Constants.u,
                "g": constants.Constants.u*1000.0,
                "m_e": constants.Constants.u/constants.Constants.m_e,
                "m_p": constants.Constants.u/constants.Constants.m_p}


# In[16]:


## print the adjaceny list 
def adjaceny(c):
    adj_list=[]

    for i in range(n):
        adj_list.append([])
    r=l_unit_dict["ang"]
    r_dict={}
    atoms=x_dict["atoms"]
    for atom in atoms:
        if atom not in r_dict:
            r_dict[atom] = constants.Constants.r_dict[atom]*r
    for i in range(n):
        for j in range(i+1,n):
            d_max = r_dict[atoms[i]] + r_dict[atoms[j]] + c*r
            dist,dist1 = distance(i,j)
            if dist < d_max:
                #x=z_dict["atoms"][i]
                #y=z_dict["atoms"][j]
                adj_list[i].append(j)
                adj_list[j].append(i)
                #adj_list[x].append(y)
                #adj_list[y].append(x)
    return adj_list
    #if atom not in r_dict_old.keys():


# In[17]:


print("atom_list:",z_dict["atoms"],"\n","connectivity_info:",adjaceny(0.4))


# In[18]:


connectivity_dict={}
a=adjaceny(0.4)


# In[19]:


c_list=[]
for i in range(n):
    temp=[]
    for j in range(len(a[i])):
        x=a[i][j]
        y=z_dict["atoms"][x]
        temp.append(y)
    c_list.append(temp)
print(c_list)


# In[20]:


for i in range(0,n):
    a=adjaceny(0.4)
    print(z_dict["atoms"][i],c_list[i])
    connectivity_dict[z_dict["atoms"][i]]=c_list[i]


# In[21]:


print_dict(connectivity_dict)


# In[22]:


## do the index matching with only connected atoms
## 


# In[23]:


#for i in range(0,n):
    #print(a[i])
#    y=[x for x in a[i]]
#    print(y)


# In[ ]:





# In[24]:


## Calculate the bond angle between three points
def angle(a,b,c):
    abc_list = [a,b,c]
    for i in range(3):
        for j in range(i+1,3):
            if abc_list[i] == abc_list[j]:
                break
    xyz = np.transpose(x_dict["xyz"])
    a_array = xyz[a]
    b_array = xyz[b]
    c_array = xyz[c]
    
    ## 3 dimensional vector p, difference of vectors a and b
    p_array = a_array - b_array
    r_array = c_array - b_array
    
    ## dot product between p and r
    p_r_dot_product = sum(p_array*r_array)
    p_norm = sum(p_array**2)**0.5
    r_norm = sum(r_array**2)**0.5
    ## The cosine of the angle between vectors p amd r in radians
    cos_t_rad = p_r_dot_product/(p_norm*r_norm)
    
    ## Angle t in radians
    t_rad = np.arccos(cos_t_rad)
    ## Converting angle t to degress
    t_deg = t_rad*180/np.pi
    return t_deg
            
    


# In[25]:


def dihedral(a,b,c,d):
    abcd_list = [a, b, c, d]
    for i in range(4):
        for j in range(i + 1, 4):
            if abcd_list[i] == abcd_list[j]:
                raise ValueError("a, b, c and d should be different.")
    abc = abs(angle(a, b, c))
    bcd = abs(angle(b, c, d))
    if not 0.0 < z < 10.0:
        raise ValueError("z should be a positive number below 10.")
    if abc < z or abc > 180.0 - z or bcd < z or bcd > 180.0 - z:
        raise ValueError("Undefined dihedral angle.")
        # The n-by-3 dimensional matrix of atomic coordinates.
    xyz = np.transpose(x_dict["xyz"])
        # 3 dimensional vectors a, b, c and d representing the atomic
        # coordinates of the corresponding atoms.
    a_array = xyz[a]
    b_array = xyz[b]
    c_array = xyz[c]
    d_array = xyz[d]
        # 3 dimensional vector qi, the difference of vectors b and a.
        # qi = b - a.
    qi = b_array - a_array
        # 3 dimensional vector qj, the difference of vectors c and b.
        # qj = c - b.
    qj = c_array - b_array
        # 3 dimensional vector qk, the difference of vectors d and c.
        # qk = d - c.
    qk = d_array - c_array
        # Vector ri, the cross product of vectors qi and qj.
    ri = np.cross(qi, qj)
        # Vector ni is vector ri divided by its norm.
    ni = ri/sum(ri**2)**0.5
        # Vector rj, the cross product of vectors qj and qk.
    rj = np.cross(qj, qk)
        # Vector ui is vector rj divided by its norm.
    ui = rj/sum(rj**2)**0.5
        # Vector uk is vector qj divided by its norm.
    uk = qj/sum(qj**2)**0.5
        # Vector uj, the cross product of vectors uk and ui.
    uj = np.cross(uk, ui)
        # x, the dot product of vectors ni and ui.
    x = sum(ni*ui)
        # y, the dot product of vectors ni and uj.
    y = sum(ni*uj)
        # The negative of the 2-argument arctangent of y/x.
    f_rad = -np.arctan2(y, x)
        # Converting angle t to degrees.
    f_deg = f_rad*180.0/np.pi
    return f_deg


# In[26]:


print("atom_list:",z_dict["atoms"],"\n")
print("connectivity_info:",adjaceny(0.4),"\n")


# In[27]:


## create the connectivity dictionary


# In[28]:


## calculate the bond angle between three points


# In[ ]:





# In[29]:


print_dict(connectivity_dict)


# In[ ]:





# In[30]:


#for i in connectivity_dict.keys():
  #  for j in range(len(connectivity_dict[i])):


# In[31]:


print(n)


# In[32]:


dist_list_temp=[]
dist_list_temp_dict={}
temp_list=list(connectivity_dict.keys())
#print(temp_list)
list1=[]
j_list=[]
k_list=[]
for i in range(len(temp_list)):
    #print(temp_list[i],connectivity_dict[temp_list[i]]
    j=int(temp_list[i].split('-')[0])
    for values in connectivity_dict[temp_list[i]]:
        k=int(values.split('-')[0])
        dist,dist1=distance(j,k)
       # print(j,k,'\t',dist,dist1)
        dist_list_temp.append(dist1)
        #print(x,":",y,'\t',j,":",k)
        l1=str(j)+":"+str(k)
        j_list.append(str(j,))
        k_list.append(str(k,))
        dist_list_temp_dict[l1]=dist1


# In[33]:


list1=list(zip((j_list),(k_list)))


# In[34]:


x=np.asarray(list1)


# In[35]:


temp_list=[]
u = np.sort(x, axis=1)

_, idx = np.unique(u, axis=0, return_index=True)

print(x[idx])
arr_list=x[idx].tolist()


# In[36]:


print(arr_list)


# In[37]:


## list comprehension of arr_list
temp_list=[]
for i in range(len(arr_list)):
    y=arr_list[i][0]+':'+arr_list[i][1]
    print(y)
    temp_list.append(y)
print(temp_list)


# In[38]:


bond_distance={}
for keys in dist_list_temp_dict.keys():
    i,j=keys.split(':')
    a= str(z_dict["atoms"][int(i)])+":"+str(z_dict["atoms"][int(j)])
    if keys in temp_list:
        bond_distance[a]=dist_list_temp_dict[keys]


# In[39]:


print(temp_list)


# In[40]:


print_dict(bond_distance)
print("\n")
print(dist_list_temp_dict)


# In[41]:


## 1. how to remove the swapped duplicates from the bond_distance list
## create the angle pair list (i,j,k)


# In[42]:


angle_list_temp=[]
angle_list_temp_dict={}
a1_list=[]
for connec in connectivity_dict.keys():
    a_list=[]
    i=int(connec.split('-')[0])
    if len(connectivity_dict[connec])!=1:
        print(connec,connectivity_dict[connec])
        a1_list.append(i)
        angle_list_temp_dict[connec]=connectivity_dict[connec]
        for values in connectivity_dict[connec]:
            j=values.split('-')[0]
            a_list.append(j)
        angle_list_temp.append(a_list)


# 

# In[43]:


print_dict(angle_list_temp_dict)


# In[ ]:





# In[44]:


## print all possible combinations of pairs
pair_list=[]
pair_list_dict={}
for i in range(len(angle_list_temp)):
    atom=[j for j in angle_list_temp_dict.keys()]
    print(atom[i])
    x=angle_list_temp[i]
    res=[[a,b] for idx, a in enumerate(x) for b in x[idx + 1:]]
    print(res)
    pair_list.append(res)
    pair_list_dict[atom[i]]=res


# In[45]:


print(pair_list_dict)


# In[ ]:





# In[46]:


l1=list(angle_list_temp_dict.keys())
print(l1)


# In[47]:


ang_pair_list=[]
ang_pair_dict={}
for i in range(len(pair_list)):
    ang_pair_temp=[]
    for j in range(len(pair_list[i])):
        x,y=pair_list[i][j][0],pair_list[i][j][1]
        z_list=str(z_dict["atoms"][int(x)])+":"+str(l1[i])+":"+str(z_dict["atoms"][int(y)])
        ang_pair_temp.append(z_list)
    ang_pair_list.append(ang_pair_temp)


# In[48]:


#ang_pair_list=[]
#ang_pair_dict={}
#for i in range(len(pair_list)):
#    for x in angle_list_temp_dict.keys():
#        print(z_dict["atoms"][int(x.split('-')[1])])
#        ang_pair_temp=[]
    #print(z_dict["atoms"][i])
#        for j in range(len(pair_list[i])):
#            x,y=pair_list[i][j][0],pair_list[i][j][1]
#            z_list=str(z_dict["atoms"][int(x)])+":"+str(z_dict["atoms"][int(x.split('-')[1])])+":"+str(z_dict["atoms"][int(y)])
#            ang_pair_temp.append(z_list)
#        ang_pair_list.append(ang_pair_temp)
        #print(pair_list[i][j], pair_list[i][j][0],pair_list[i][j][1],z_dict["atoms"])
                                                        


# In[49]:


print(ang_pair_list)


# In[50]:


print(ang_pair_list)


# In[51]:


## Make the angle pair combinations 
angle_list_temp1=[]
for i in range(len(angle_list_temp)):
    temp1=[]
    #print(i,pair_list[i])
    for j in range(len(pair_list[i])):
                   print(pair_list[i][j][0],l1[i].split('-')[0],int(pair_list[i][j][1]))
                   angle_temp=angle(int(pair_list[i][j][0]),int(l1[i].split('-')[0]),int(pair_list[i][j][1]))
                   print(angle_temp)
                   temp1.append(angle_temp)
    angle_list_temp1.append(temp1)


# In[52]:


print(angle_list_temp1)


# In[53]:


## create the keys for angle_list_dict
angle_dict1={}
ang_l1=[]
for i in range(len(ang_pair_list)):
    #print(ang_pair_list[i])
    #print(angle_list_temp1[i])
    for j in range(len(ang_pair_list[i])):
        print(ang_pair_list[i][j],angle_list_temp1[i][j])
        angle_dict1[str(ang_pair_list[i][j])]=angle_list_temp1[i][j]


# In[54]:


print(connectivity_dict)


# In[55]:


print(len(angle_dict1))


# In[56]:


## check for nitrogen
sum=0
for atom in angle_dict1.keys():
    i=atom.split(":")
    #print(i,i[1].split('-')[0])
    if (i[1].split('-')[1])=='N':
        sum=sum+angle_dict1[atom]
        #print(atom, angle_dict1[atom])


# In[57]:


print("Connectivity")
print_dict(connectivity_dict)


# In[58]:


print("Bond Distance")
print_dict(bond_distance)


# In[59]:


print("Angles")
print_dict(angle_dict1)
print('Sum of angles around N',sum)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


from pathlib import Path


# In[3]:


file_path = Path(r"Propane-new.xyz")


# In[ ]:





# In[ ]:





# In[4]:


def open_xyz(path, encoding="utf-8"):
    with open(path, encoding=encoding) as file:
        read= file.read().strip()
       # read=file.readlines()
        #return [read,5]
        rows=read.split("\n")
        data=[rows[0],rows[1],rows[2]]
        data=data+[row for row in rows[2:] if row.strip !=""]
        return rows,data
    


# In[5]:


#file=open_xyz("Propane-new.xyz")
rows,data=open_xyz("Propane-new.xyz")


# In[ ]:





# In[6]:


def open_log(path, encoding="utf-8"):
    with open(path,encoding=encoding) as file:
        lines=file.readlines()
        return lines


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





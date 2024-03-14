##################################################################
#Dodelavat komentare
#Cteni dat z exif z obrazku
#Filtrovani krivek (namisto resample udelat prumerovani na nejakem poctu vzorku)
#Integrace krivky (udela scipy sims)
#Vycisteni krivek od stejnych hodnot v jedne ose
#Integrace bambuli (Hanka/Ja)
##################################################################

##########################################################################################
#Prepare database and read data from txt files
##########################################################################################

def read_keys_from_dir_KROUPA(id_dir,extensions):
    #read names of files with extensions into list
    #input
    #id_dir     = string, where can be files found
    #extensions = list of strings, extensions with dot of files to use ['.TRA','.DAT']
    #output
    #list_keys  = list of strings with names of files with name of subdirectories if they were in any and no extension of the file
    
    import os

    ok_code   = bool(1)
    
    list_keys = []

    if os.path.isdir(id_dir):

        for path,dirs,files in os.walk(id_dir):
    
            for id_file in files:
            
                use_ext = bool(0)
            
                for ext in extensions:
                
                    if ext in id_file:
                    
                        use_ext = bool(1) 
            
                if use_ext:
                
                    key_total = os.path.join(path, id_file)
                    
                    key       = key_total.replace(id_dir,'') 
                    
                    list_keys = list_keys + [key[:-len(ext)]]
                    
    else:                    

        print '      *read_keys_from_dir_KROUPA - directiroy = ',id_dir,' does not exist!'
        ok_code = bool(0)

    return list_keys,ok_code 

def create_database_from_keys_KROUPA(keys):
    #from list of strings create dictionary with kyes from list
    #input
    #keys     = list of keys
    #output
    #d        = ordereddict
    #ok_code  = bool, if false anything is in list
    
    import collections
    
    d = collections.OrderedDict()
    
    ok_code = bool(1)
    
    if len(keys)<1:
    
        print '      *create_database_from_keys_KROUPA - list of keys is emepty'
        ok_code = bool(0)  

    else:
    
        for key in keys:
        
            d[key]  = collections.OrderedDict() 

    return d,ok_code
    
def read_TXTs_KROUPA(dir,keys,exts,no_lines_to_skip,data_types,columns):
    #read data from text files with given lines to skip and from given columns from directory with names as keys
    #input
    #dir    = directory where are the text files, it should end with "\\"
    #keys   = list of strings, names of text files
    #exts   = where id extension of given file (it is necessary, if you have two files for one measurement)
    #no_lins_to_skip  = dictionary of lists of one integer, number lines to skip in each experiement
    #data_types = data type of data in columns (usualy float or str) - what i supported by loadtxt can be used
    #columns          = dictionary of lists of one integer, from which column read the data 
    #output
    #data     = dictionary with only read data
    #ok_code  = 
    
    import collections
    import numpy as np
    import os

    ok_code = bool(1)
    
    data  = collections.OrderedDict()
    
    for key in keys:
        
        id_file = dir+key+exts[key][0]
        if  os.path.exists(id_file):
        
            print '      *read_TXTs_KROUPA - reading ',id_file
            if (key in data_types.keys()) and (key in columns.keys()):
                f = open(id_file, 'r')
                data[key]  = np.loadtxt(f, dtype=data_types[key][0], comments='#', delimiter=None, converters=None, skiprows=no_lines_to_skip[key][0], usecols=columns[key], unpack=False, ndmin=0)
                f.close()
            else:
                print '      *read_TXTs_KROUPA - key "',key,'" does not exists in data_types or columns'
                ok_code = bool(0)            
            
        else:
            
            print '      *read_TXTs_KROUPA - file "',id_file,'" does not exists '    
            ok_code = bool(0) 
            
    return data,ok_code
    
def find_indices_for_column_TXT_KROUPA(quantity_TRA,field):
    #find index of column in data files if name of field is given and is in field quantity_TRA
    #input
    #quantity_TRA = list of strings, where in the dictionary can be found specifications of columns
    #field        = list of one string, name of column in quantity_TRA
    #output
    #indices      = dictionary, key is name of experiment and data are list with integer defining number of column
    #ok_code      = boolean

    import collections

    ok_code = bool(1)
    
    indices = collections.OrderedDict()
    
    for key_quantity in quantity_TRA.keys():
    
        field_found = bool(0)
        for i in range(0,len(quantity_TRA[key_quantity])):
            if quantity_TRA[key_quantity][i]==field:
                index = i
                field_found = bool(1)
      
        if field_found:
            indices[key_quantity] = [index]
        else:
            print '      *find_indices_for_column_TRA_KROUPA - key "',key_quantity,'" - field "',field,'" was not found in "',quantity_TRA[key_quantity],'"'    
            ok_code = bool(0)            

    return indices,ok_code
    
##########################################################################################
#Extract data and find unique data in fields
##########################################################################################

def find_unique_data_KROUPA(d,keys,field):
    #finds unique values in given filed in d (works only fo one value i suppose)
    #input
    #d      = ordereddictionary,data
    #keys   = list of strings,keys in which search
    #field  = name of field where to search for unique data
    #ouput
    #unique_data  = list of unique data
    #ok_code      = false if no data found 

    ok_code = bool(1)

    possible_data,ok_code = extract_data_fields_KROUPA(d,keys,field)

    unique_data = []

    if len(possible_data)>0:

        for key in possible_data.keys():
        
            if not (possible_data[key][0] in unique_data):
            
                unique_data = unique_data + possible_data[key]
                
    else:
    
        print '      *find_unique_data_KROUPA - in field =',field,' are no data'
        ok_code = bool(0)                 

    return  unique_data,ok_code

def extract_data_field_KROUPA(d,field):
    #extract data from given field
    #input
    #d          = ordereddictioanry
    #key path   = list of strings, name of field from where to extract data ['tensile','displacement']
    #output
    #data       = everything what is in the field
    #ok_code    = false if 

    ok_code = bool(1)

    for key in field:

        d = d[key]

    return d,ok_code

def extract_data_fields_KROUPA(d,keys,field):    

    import collections
    
    ok_code = bool(1) 
    
    data  = collections.OrderedDict() 

    for key in keys:
        
        if check_existence_field_KROUPA(d[key],field):
        
            data[key],ok_code_field   = extract_data_field_KROUPA(d[key],field)
            
        else:
        
            print '      *extract_data_fields_KROUPA - Field = ',field,', Data for key = ',key,' do not exist!'
            ok_code = bool(0)
            
    if len(data.keys())<1:
    
        ok_code = bool(0)                        
        
    return data,ok_code
    
def extract_number_in_vector_fields_KROUPA(d,indices):

    import collections

    ok_code = bool(1)
    
    numbers = collections.OrderedDict()
    
    for key in d.keys():
    
        if key in indices.keys():
            numbers[key]  = [d[key][indices[key][0]]]
        else:                    
            print '      *extract_number_in_vector_fields_KROUPA - key = ',key,', not exist in indices'
            ok_code = bool(0)

    return numbers,ok_code    

def arithmetic_data_fields_KROUPA(d,numbers,operation):    

    import numpy as np
    
    ok_code = bool(1) 

    for key in d.keys():
        
        if operation[0]=='+':
        
            d[key]   =  np.array(d[key])  + numbers[key][0]
            
        elif operation[0]=='-':
        
            d[key]   =  np.array(d[key])  - numbers[key][0]
        
        elif operation[0]=='*':

            d[key]   =  np.array(d[key])  * numbers[key][0]          

        elif operation[0]=='/':
        
            d[key]   =  np.array(d[key])  / numbers[key][0]         
        
    return d,ok_code
    
def fill_data_default_KROUPA(d,keys,default):    
    #fill field with default data
    #input
    #d        = ordereddict,database
    #keys     = list of keys for check and/or fill   
    #default  = what to add
    #output
    #d        = database
    #ok_code  = if at least one field were fiiled ok_code=true otherwise false 

    ok_code = bool(0) 

    for key in keys:
        
        if not (key in d.keys()):
        
            ok_code   = bool(1)
            d[key]   =  default
        
    return d,ok_code
    
##########################################################################################
#Create, delete, rename keys of fileds
##########################################################################################
    
def create_dict_field_KROUPA(d,key):

    import collections
    
    ok_code  = bool(1)
    
    if not check_existence_field_KROUPA(d,key):

        d[key[0]]   = collections.OrderedDict()

    else:
    
        print '      *create_dict_field_KROUPA - Key ',key,' exist!'
        ok_code  = bool(0)                
    
    return d,ok_code 
    
def rename_field_KROUPA(d,key_old,key_new):

    ok_code  = bool(1)
    
    if check_existence_field_KROUPA(d,key_old):

        zaloha        = d[key_old[0]]
        
        del d[key_old[0]]
        
        d[key_new[0]] = zaloha 

    else:
    
        print '      *rename_field_KROUPA - Key ',key_old,' does not exist!'
        ok_code  = bool(0)                
    
    return d,ok_code

def delete_field_KROUPA(d,key):

    ok_code  = bool(1)
    
    if check_existence_field_KROUPA(d,key):

        del d[key[0]]
        
    else:
    
        print '      *delete_field_KROUPA - Key ',key,' do not exist!'
        ok_code  = bool(0)                
    
    return d,ok_code
    
def delete_fields_KROUPA(d,keys):

    ok_code  = bool(1)
    
    for key in keys:
    
        d,ok_code_delete = delete_field_KROUPA(d,[key])
        
        if not ok_code_delete:
        
            ok_code = bool(0)

    return d,ok_code               

##########################################################################################
#Keys operations
##########################################################################################

def merge_lists_KROUPA(list1,list2):
    #merge lists
    #input
    #list1,list2    = lists, lists for intersection
    #output
    #merged_list = merged list
    #ok_code        = bool, false if intersect_list is empty
    
    ok_code     = bool(1)

    merged_list = list(set(list1 + list2))
    
    if len(merged_list)<1:
    
        ok_code = bool(0)

    return merged_list,ok_code        

def intersect_lists_KROUPA(list1,list2):
    #intersect lists
    #input
    #list1,list2    = lists, lists for intersection
    #output
    #intersect_list = intersected list
    #ok_code        = bool, false if intersect_list is empty
    
    ok_code         = bool(1)
    
    intersect_list  = list(set(list1) & set(list2))

    if len(intersect_list)<1:
    
        ok_code = bool(0)
    
    return intersect_list,ok_code
    
def list_not_in_list_KROUPA(list1,list2):
    
    ok_code         = bool(1)
    
    filtered_list = []
    
    for key in list2:
    
        if not (key in list1):
        
            filtered_list = filtered_list + [key]
            
    if len(filtered_list)<1:
    
        ok_code = bool(0)                
    
    return filtered_list,ok_code

def filter_keys_key_KROUPA(keys,condition_string,condition_type):
    #return list of keys which satisfied a condition
    #keys = list of strings (keys), ['prdlavka.TRA','tahovka_tkanicky.TRA','casovy_zaznam_tahovky.DAT']
    #condition_string = list of one string, string which is tested for existence
    #condition_type   = list of one special character, ['=']-string is equal or ['in']-string is in name of key
    #output
    #keys_filtered    = list of strings,list of keys
    #ok_code          = bool,indicates wrong condition type or no data in keys_filtered
    
    ok_code = bool(1)

    keys_filtered = []

    for key in keys:
        
        use_key = bool(0)
        
        if condition_type[0]=='=':
        
            for cs in condition_string:  
                
                if cs==key:
            
                    use_key = bool(1)              

        elif condition_type[0]=='in':
        
            for cs in condition_string:
                
                if cs in key:
            
                    use_key = bool(1)
        
        else:
        
          print '      *filter_keys_key_KROUPA - you are using not supported condition_type - supported types: [''=''],[''in'']'
          ok_code = bool(0)
          
        if use_key:
        
            keys_filtered = keys_filtered + [key]
            
    if len(keys_filtered)<1:
    
        print '      *filter_keys_key_KROUPA - no data in output!'
        ok_code = bool(0)            
        
    return keys_filtered,ok_code
    
def filter_keys_field_KROUPA(d,field):
    #Filters keys by existence of field
    
    ok_code = bool(1)

    keys_filtered = []

    for key in d.keys():
        
        use_key = bool(0)
        
        if check_existence_field_KROUPA(d[key],field):
                
            use_key = bool(1)
        
        else:
        
          print '      *filter_keys_field_KROUPA - field ',field,' does not exist in key: ',key
          ok_code = bool(0)
          
        if use_key:
        
            keys_filtered = keys_filtered + [key]
            
    if len(keys_filtered)<1:
    
        ok_code = bool(0)            
        
    return keys_filtered,ok_code    
    
def filter_keys_data_KROUPA(d,condition_field,condition_data,condition_type):
    #Filters keys by condition for the data in field
    #First it checks existence of key
    #Input
    #condition_field = list, name of field which is used for filtering of keys
    #condition data = list, what data are udsed for filtering   
    #condition type = what condition is used possible conditions are =, >,<, in
    #output
    #keys_filtered = list of kyes which satisfied condition

    ok_code = bool(1)

    keys_filtered = []

    for key in d.keys():
        
        if not check_existence_field_KROUPA(d[key],condition_field):
        
            print '      *filter_keys_data_KROUPA - condition_field: ',condition_field,' does not exist in field with key: ',key
            ok_code = bool(0)     
            
        else:
        
            data,ok_code_extract  = extract_data_field_KROUPA(d[key],condition_field)
        
            use_key = bool(0)
            
            if condition_type[0]=='=':
                
                if data==condition_data:
                
                    use_key = bool(1)              

            elif condition_type[0]=='<':
            
                if data[0]<condition_data[0]:
                
                    use_key = bool(1)

            elif condition_type[0]=='>':
            
                if data[0]>condition_data[0]:
                
                    use_key = bool(1)

            elif condition_type[0]=='in':
            
                if condition_data[0] in data[0]:
                
                    use_key = bool(1)
            
            else:
            
              print '      *filter_keys_data_KROUPA - you are using not supported condition_type - supported types: [''='']-strings and numbers also lists;[''>''],[''<'']-numbers - only one number should be used'
              ok_code = bool(0)
              
            if use_key:
            
                keys_filtered = keys_filtered + [key]
        
    return keys_filtered,ok_code

def print_dict_tree_KROUPA(d,level):
    #recursively print data tree
    
    empty_string=''
    for i in range(0,7*level):
        empty_string=empty_string+' '
    for key in d.keys():
        print "%s[%s]" % (empty_string,key)
        if isinstance(d[key], dict):
            print_dict_tree_KROUPA(d[key],level+1)
    
def check_existence_field_KROUPA(d,key_path):
    #check existance of field
    #input
    #d        = ordereddict,database
    #key_path = list of key, key path, ['tensile','displacement']
    #output
    #ok_code  = bool, true if field exists

    ok_code = bool(0)
    
    for key in d.keys():
        
        if key_path[0]==key:
        
            ok_code = bool(1)
            
    if ok_code:
        length_kyes_check = len(key_path)
        if length_kyes_check>1:
            ok_code = check_existence_field_KROUPA(d[key_path[0]],key_path[1:length_kyes_check])
            
    return ok_code

##########################################################################################
#working and mathematical functions
##########################################################################################

def extract_column_list_of_lists_KROUPA(data_RAW,column):
    #extract column from list of lists array
    #input
    #data_RAW = lists of lists
    #column   = integer,index of column to extract
    #output
    #data     = numpy array, extracted vector
    #ok_code  = bool, false if bad column index is given

    ok_code = bool(1)

    c = len(data_RAW[1])

    if not column>=c:

        data  = [item[column] for item in data_RAW]
        
    else:
    
        print '      *extract_column_list_of_lists_KROUPA - column number ',column,'is greater than number of columns in data_RAW ',c 
        ok_code=bool(0)
    
    return data,ok_code

def extract_column_array_KROUPA(data_RAW,column):
    #extract column from numpy array
    #input
    #data_RAW = numpy aray
    #column   = integer,index of column to extract
    #output
    #data     = numpy array, extracted vector
    #ok_code  = bool, false if bad column index is given
    
    import numpy as np

    ok_code = bool(1)

    c = np.size(data_RAW,1)

    if not column>=c:

        data  = data_RAW[:,column]
        
    else:
    
        print '      *extract_column_array_KROUPA - column number ',column,'is greater than number of columns in data_RAW ',c 
        ok_code=bool(0)
    
    return data,ok_code
    
def extract_raw_list_of_lists_KROUPA(data_RAW,raw):
    #extract raw from list opf lists
    #input
    #data_RAW = list of lists
    #column   = integer,index of raw to extract
    #output
    #data     = numpy array, extracted vector
    #ok_code  = bool, false if bad raw index is given

    ok_code = bool(1)

    r = len(data_RAW[0])

    if not raw>=r:

        data  = data_RAW[raw]
        
    else:
    
        print '      *extract_raw_list_of_lists_KROUPA - raw number ',raw,'is greater than number of raws in data_RAW ',r 
        ok_code=bool(0)
    
    return data,ok_code     
    
def extract_raw_array_KROUPA(data_RAW,raw):
    #extract raw from numpy array
    #input
    #data_RAW = numpy aray
    #column   = integer,index of raw to extract
    #output
    #data     = numpy array, extracted vector
    #ok_code  = bool, false if bad raw index is given
    
    import numpy as np

    ok_code = bool(1)

    r = np.size(data_RAW,0)

    if not raw>=r:

        data  = data_RAW[raw,:]
        
    else:
    
        print '      *extract_raw_array_KROUPA - raw number ',raw,'is greater than number of raws in data_RAW ',r 
        ok_code=bool(0)
    
    return data,ok_code     
    
def find_max_point_KROUPA(x,y):    

    import numpy as np
    
    ok_code     = bool(1)
    
    if len(y)>0:
    
        max_y       =   np.array([max(y)])

        index_found = bool(0)
        index       = 0
        while not index_found:
        
            if max_y>y[index]:
                
                index = index+1
                
            else:
                
                index_found = bool(1)
                max_x       = np.array([x[index]])                
        
    else:
    
        max_y       =   np.array([])        
        max_x       =   np.array([])
        index       =   0
        ok_code     = bool(0)    
    
    return max_x,max_y,index,ok_code

def find_interval_KROUPA(x,y,bounds):    
    #finds data x_i and y_i in x and y which x data are in interval given by bounds
    #input 
    #x,y    = numpy array or lists of data
    #bounds = numpy array or list with lower and upper bounds, example = [0,1]
    #output
    #x_i,y_i  = numpy array or lists of data, data in interval bounds
    #ok_code  = false if there are no data in interval
    
    import numpy as np

    ok_code = bool(1)

    x_b  = x[x>=bounds[0]]
    y_b  = y[x>=bounds[0]]
    
    x_i  = x_b[x_b<=bounds[1]]
    y_i  = y_b[x_b<=bounds[1]]
    
    if len(x_i)<1 or len(y_i)<1:
    
        print '      *find_interval_KROUPA - there are no data in interval'
        ok_code = bool(0) 
        
    return x_i,y_i,ok_code
    
def find_neighbors_KROUPA(x,y,point):
    #find neibors for value point in x dimension and returns x and y coordinates of neibors (left,right)
    #input
    #x      = numpy array,list,vector
    #y      = numpy array,list,vector
    #point  = number, float, point in x dimension around which is neibors searched
    #output
    #np.array([b_x[-1],t_x[0]]) = vector in x dimension of two neibors
    #np.array([b_y[-1],t_y[0]]) = vector in y dimension of two neibors
    #ok_code  = bool, is flase if no neibors are found
    
    import numpy as np
    
    ok_code = bool(1)
    
    b_x = x[x<=point]
    b_y = y[x<=point]
    
    t_x = x[x>=point]
    t_y = y[x>=point]

    if  b_x[-1]==t_x[0]:

        if  point<=min(x):
        
            t_x  = x[x>point]
            t_y  = y[x>point]
            
        else:  
        
            b_x  = x[x<point]
            b_y  = y[x<point]        
    
    if len(t_x)<1 or len(t_y)<1 or len(b_x)<1 or len(b_y)<1:
    
        ok_code = bool(0)
    
    return np.array([b_x[-1],t_x[0]]),np.array([b_y[-1],t_y[0]]),ok_code
    
def resample_curves_KROUPA(data_x,data_y,points):
    #resample curves
    #input
    #data_x,data_y  = vectors of data
    #points         = vector with x coord. where to resample
    #output
    #data_resampled = resampled data in y directions (recalculated y values for points values)
    #ok_code        = false if no data for averaging exists
    
    import collections
    import numpy as np
    
    ok_code = bool(1)         
    
    if len(data_x)!=len(data_y):
        print '      *resample_curves_KROUPA - data are not the same length'
        ok_code = bool(0)
        
    if min(data_x)>points[0]:
        print '      *resample_curves_KROUPA - points are too large (data_x[0]>points[0])'
        ok_code = bool(0)

    if max(data_x)<points[-1]:
        print '      *resample_curves_KROUPA - points are too large (data_x[-1]<points[-1])'
        ok_code = bool(0)
        
    data_resampled  = np.array([])

    if ok_code:

        for i in range(0,len(points)):
    
            nb_x,nb_y,ok_code_find      = find_neighbors_KROUPA(data_x,data_y,points[i])
            
            polynom                     = np.polyfit(nb_x,nb_y,1)
            
            data_resampled  = np.append(data_resampled,np.polyval(polynom,points[i]))
        
    return data_resampled,ok_code         

def find_indices_local_extremes_KROUPA(x, delta):
    #find indices of local extremes on single vectors - it is assumed that minimum is between two maximas and maximum is between two minmas
    #input
    #x      = numpy array or list,vector,
    #delta  = list of one numebr,parameter defines weather is the point local extreme. If difference of two values of x is greater than delta and index is after the extreme, than we find extreme
    #output
    #ind_min,ind_max  = lists of integers,indices of extremes  
    
    ind_max = []
    ind_min = []
      
    if (len(delta))>1 or (len(delta))<1:
      print '      *find_indices_local_extremes_KROUPA - Input argument delta ',str(delta),' must be a scalar'
    
    if delta<=0:
      print '      *find_indices_local_extremes_KROUPA - Input argument delta ',str(delta),' must be positive';
    
    mn    =  1e64 
    mx    = -1e64
    
    lookformax = bool(1)
    
    for i in range(0,len(x)):

      this  = x[i]

      if this > mx:
          mx      = this
          ind_mx  = i
      if this < mn:
          mn = this
          ind_mn  = i
      
      if lookformax:
        if this < mx-delta:
          ind_max = ind_max + [ind_mx]
          mn = this
          lookformax = bool(0)
      else:
        if this > mn+delta:
          ind_min = ind_min + [ind_mn]
          mx = this
          lookformax = bool(1)

    return ind_max, ind_min
    
def find_int_SRBOVA(x,y,ind_min,ind_max):
    #Finds intersections of bambule (two curves with one intersection) between index of maximum a index of maximum+1 - there mus be index of minimum between them
    #In order to work well - zero index is added to begin of the ind_min vector 
    #input
    #x  = vector of x data 
    #y  = vector of y data
    #ind_min  = vector of indices of local minimums
    #ind_max  = vector of indices of local maximums
    #output
    #x_int, y_int = vektors of coordinates of intersection x and y 
    
    import numpy as np
    
    ok_code   = bool(1)
    
    ind_min = [0] + ind_min
    
    ind_max = ind_max + [len(y)-1]

    x_int =  np.array([])
    y_int =  np.array([])
    
    min_min = int(min([len(ind_max),len(ind_min)]))
     
    for ind in xrange(min_min-1):
    
        print '            Intersection = ',ind+1,'/',min_min-1,' done'
        
        if (ind_max[ind] > ind_min[ind+1]) or (y[ind_max[ind]] < y[ind_min[ind+1]]) or (x[ind_max[ind]] < x[ind_min[ind+1]]):
            ok_code = bool(0)
            print '      *find_int_SRBOVA - chyba v klesajici casti cyklu'
            if (ind_max[ind] > ind_min[ind+1]):
                print '            *chyba v klesajici casti cyklu 1: ind_max[ind] > ind_min[ind+1]: '
                print '                  ind_max[ind] ', ind_max[ind]
                print '                  ind_min[ind+1] ', ind_min[ind+1]
            if (y[ind_max[ind]] < y[ind_min[ind+1]]):
                print '            *chyba v klesajici casti cyklu 2: y[ind_max[ind]] < y[ind_min[ind+1]]: '
                print '                  y[ind_max[ind]] ', y[ind_max[ind]]
                print '                  y[ind_min[ind+1]] ', y[ind_min[ind+1]]
            if (x[ind_max[ind]] < x[ind_min[ind+1]]):
                print '            *chyba v klesajici casti cyklu 3: x[ind_max[ind]] < x[ind_min[ind+1]]: '
                print '                  x[ind_max[ind]] ', x[ind_max[ind]]
                print '                  x[ind_min[ind+1]] ', x[ind_min[ind+1]]
           
        if (ind_min[ind+1] > ind_max[ind+1]) or (y[ind_min[ind+1]] > y[ind_max[ind+1]]) or (x[ind_min[ind+1]] > x[ind_max[ind+1]]):
            ok_code = bool(0)
            print '      *find_int_SRBOVA - chyba v rostouci casti cyklu'
            if (ind_min[ind+1] > ind_max[ind+1]):
                print '            *chyba v rostouci casti cyklu 1: ind_min[ind+1] > ind_max[ind+1]:'
                print '                  ind_min[ind+1] ', ind_min[ind+1]
                print '                  ind_max[ind+1] ', ind_max[ind+1]
            if (y[ind_min[ind+1]] > y[ind_max[ind+1]]):
                print '            *chyba v rostouci casti cyklu 2: y[ind_min[ind+1]] > y[ind_max[ind+1]]:'
                print '                   y[ind_min[ind+1]] ', y[ind_min[ind+1]]
                print '                   y[ind_max[ind+1]] ', y[ind_max[ind+1]]
            if (x[ind_min[ind+1]] > x[ind_max[ind+1]]):
                print '             *chyba v rostouci casti cyklu 3: x[ind_min[ind+1]] > x[ind_max[ind+1]]:'
                print '                   x[ind_min[ind+1]] ', x[ind_min[ind+1]]
                print '                   x[ind_max[ind+1]] ', x[ind_max[ind+1]]
    
        found_int = 0
        
        cast_pr_int = int(0.1 * len(x[ind_min[ind]:ind_max[ind]]))
        
        x_down  = x[ind_max[ind] - cast_pr_int  : ind_min[ind+1] + 1]
        y_down  = y[ind_max[ind] - cast_pr_int  : ind_min[ind+1] + 1]
        x_up    = x[ind_min[ind+1]              : ind_max[ind+1]+ 1 + int(0.1*len(x_down))]
        y_up    = y[ind_min[ind+1]              : ind_max[ind+1]+ 1 + int(0.1*len(x_down))]
        
        for j in xrange(len(x_down)-1):
            p0 = np.array([x_down[j], y_down[j]])
            p1 = np.array([x_down[j+1], y_down[j+1]])
            for k in xrange(len(x_up)-1):
                            
                q0 = np.array([x_up[k], y_up[k]])
                q1 = np.array([x_up[k+1], y_up[k+1]])
                                        
                try:
                     params = np.linalg.linalg.solve(np.column_stack((p1-p0, q0-q1)), q0-p0)                     
                except np.linalg.linalg.LinAlgError as err:
                    if 'Singular matrix' in err.message:
                        continue
                if np.all((params >= 0) & (params <= 1)):
                    inter = p0 + params[0]*(p1 - p0)
                    x_int = np.append(x_int,inter[0])
                    y_int = np.append(y_int,inter[1])
                    found_int = 1
                    
                if bool(found_int):
                    break
            if bool(found_int):
                break
                
    return x_int, y_int, ok_code 
    
def integrate_cycles_area_KROUPA(cycles):
    #integrate area in cycles of one curve
    #input
    #cycles = dictionary with cycles
    #output
    #cycles_area = dictionary with cycles area

    import numpy as np
    import collections
    from scipy.integrate import simps

    ok_code = bool(1)

    cycles_area = collections.OrderedDict()

    cycles_area['A_1']  = np.zeros(len(cycles.keys()))
    cycles_area['A_2']  = np.zeros(len(cycles.keys()))
    cycles_area['A']    = np.zeros(len(cycles.keys()))     
    
    for no_cycle in cycles.keys():
    
        cycles_area['A_1'][no_cycle]  = simps(cycles[no_cycle]['y_1'],cycles[no_cycle]['x_1'])
        cycles_area['A_2'][no_cycle]  = simps(cycles[no_cycle]['y_2'],cycles[no_cycle]['x_2'])
        cycles_area['A'][no_cycle]    = cycles_area['A_1'][no_cycle] + cycles_area['A_2'][no_cycle]      

    return cycles_area,ok_code
    
def find_cycles_KROUPA(x,y,x_int,y_int,ind_min):
    #Finds cycles and inserts each cycle into separate subtree of dictioanry
    #In order to work well - zero index is added to begin of the ind_min vector 
    #input
    #x = vector of x data
    #y = vector of y data
    #x_int, y_int = vector of coordinates of each intersection
    #ind_min = vector of indices of minimums of cycles
    #output
    #cycles = dictionary of cycles with numpy arrays of x and y data
    
    import numpy as np
    import collections
    
    ok_code =  bool(1)

    if not len(x)>0:
        print '      *find_cycles_KROUPA - len(x) = 0'
        ok_code =  bool(0)

    if not len(y)>0:
        print '      *find_cycles_KROUPA - len(y) = 0'
        ok_code =  bool(0)
        
    if not len(x_int)>0:
        print '      *find_cycles_KROUPA - len(x_int) = 0'
        ok_code =  bool(0)

    if not len(y_int)>0:
        print '      *find_cycles_KROUPA - len(y_int) = 0'
        ok_code =  bool(0)        

    if not len(ind_min)>0:
        print '      *find_cycles_KROUPA - len(ind_min) = 0'
        ok_code =  bool(0)
    
    if not len(x_int)==len(y_int):
        print '      *find_cycles_KROUPA - intersections are not the same lengths in x and y directions - len(x_int),len(y_int) = ',len(x_int),' = ',len(y_int)
        ok_code =  bool(0)    

    no_of_cycles  = min(len(x_int),len(y_int),len(ind_min))

    if no_of_cycles<1:
        print '      *find_cycles_KROUPA - found no valid cycle' 
        ok_code =  bool(0)
        
    cycles = collections.OrderedDict()

    for no_cycle in range(0,no_of_cycles):
    
        cycles[no_cycle]  =   collections.OrderedDict()
        
        cycles[no_cycle]['x_1']  =   np.array([])
        cycles[no_cycle]['y_1']  =   np.array([])                
        cycles[no_cycle]['x_2']  =   np.array([])
        cycles[no_cycle]['y_2']  =   np.array([])
    
    if ok_code:
    
        for no_cycle in range(0,no_of_cycles):

            x_1 = [x[ind_min[no_cycle]]]
            y_1 = [y[ind_min[no_cycle]]]

            x_2 = [x[ind_min[no_cycle]]]
            y_2 = [y[ind_min[no_cycle]]]
        
            #x_1 side
            ind = ind_min[no_cycle]-1
            while (x[ind]<= x_int[no_cycle]) and (y[ind]<=y_int[no_cycle]) and ind>=0:
                                        
                if not (x[ind] in x_1):
                    
                    x_1 = [x[ind]]+x_1
                    y_1 = [y[ind]]+y_1
                    
                ind = ind-1

                if ind<0:
                  break                    
                        
            #x_2 side
            ind = ind_min[no_cycle]+1
            while (x[ind]<=x_int[no_cycle]) and (y[ind]<=y_int[no_cycle]) and ind<len(x):
                                        
                if not (x[ind] in x_2):
                    
                    x_2 = x_2+[x[ind]]
                    y_2 = y_2+[y[ind]]                                       
    
                ind = ind+1

                if ind>len(x):
                      break
                
            cycles[no_cycle]['x_1']  =   np.array(x_1)
            cycles[no_cycle]['y_1']  =   np.array(y_1)                
            cycles[no_cycle]['x_2']  =   np.array(x_2)
            cycles[no_cycle]['y_2']  =   np.array(y_2)
            
    return cycles,ok_code      

#######################################################
#Performing operations on data dict()
#######################################################

def find_curves_max_min_KROUPA(d):

    import collections
    
    ok_code = bool(1)
    
    d_min = collections.OrderedDict()
    d_max = collections.OrderedDict()

    for key in d.keys():
    
        if  len(d[key])>0:
    
            d_min[key] = [min(d[key])]
            d_max[key] = [max(d[key])]
            
        else:
        
            print '      *find_curves_max_min_KROUPA - wierd data for key = ',key
            ok_code = bool(0)
        
    return  d_min,d_max,ok_code
    
def average_curves_KROUPA(data):
    #average curves stored in dictionary
    #input
    #data = ordereddictionary
    #output
    #data_averaged  = nupmy array with averaged data
    #ok_code  = false if some key is mising or if no data for averaging are present
    
    import collections
    import numpy as np
    
    ok_code = bool(1)
    
    keys = data.keys()
    
    if len(data)>0:
    
        if len(data[keys[0]])>0:
    
            max_len = 0
            for i in range(0,len(keys)):
            
                if max_len<len(data[keys[i]]):
            
                    max_len = len(data[keys[i]])    

            data_averaged     = np.zeros(max_len)
            data_no_averaged  = np.zeros(max_len) 
            
            for i in range(0,len(keys)):
                                                  
                for j in range(0,len(data[keys[i]])):
                
                    data_averaged[j]    = data_averaged[j] + data[keys[i]][j]
                    data_no_averaged[j] = data_no_averaged[j] + 1.0
            
            data_averaged = data_averaged/data_no_averaged
            
        else:
        
            print '      *average_curves_KROUPA - data for key',keys[0],' is missing!'                        
            data_averaged = np.array([])
            ok_code = bool(0)

    else:
    
        print '      *average_curves_KROUPA - no data for averaging!'                        
        data_averaged = np.array([])
        ok_code = bool(0)
                
    return data_averaged,ok_code
    
def calculate_common_base_KROUPA(data_min,data_max):
    #finds common base of curves from mminimums and maximums of the signals
    #input
    #data_min = ordereddictionary, data with minimums of signals
    #data_max = ordereddictionary, data with maximums of signals
    #output
    #min_x,max_x  = minimum and maximum of signals
    #ok_code      = is flase if no data or min_x is greater than max_x
    
    import numpy as np

    ok_code = bool(1)

    minimums  = np.array([])
    maximums  = np.array([])
    
    for key in data_min.keys():
    
        minimums  = np.append(minimums,data_min[key])
        maximums  = np.append(maximums,data_max[key])
            
    if len(minimums)>0 and len(maximums)>0:
        
        min_x = np.array([max(minimums)])
        max_x = np.array([min(maximums)])
        
        if min_x[0]>max_x[0]:
        
            print '      *calculate_common_base_KROUPA - min_x is greater that max_x for averaging!'
            ok_code = bool(0)
            
    else:
        
        print '      *calculate_common_base_KROUPA - no data for averaging!'
        min_x = np.array([])
        max_x = np.array([])
        ok_code = bool(0)
        
    return min_x,max_x,ok_code
    
def mean_std_max_min_KROUPA(d):
    #for value in dictionary calculates mean,std,max and min
    
    import numpy as np
    
    ok_code = bool(1)
    
    data  = np.array([]) 
    
    for key in d.keys():

        if len(d[key])>1 or len(d[key])<1:    
            
            print '      *mean_std_max_min_KROUPA - key =',key,' - bad data ',d[key]
            ok_code = bool(0)
            
        else:            
            
            data  = np.append(data,d[key][0])
    
    if len(data)>0:        
    
        mean_value  = [np.mean(data)]
        std_value   = [np.std(data)]
        max_value   = [np.max(data)]
        min_value   = [np.min(data)]
        
    else:
    
        print '      *mean_std_max_min_KROUPA - no data were used'
        mean_value  = [np.array([])]
        std_value   = [np.array([])]
        max_value   = [np.array([])]
        min_value   = [np.array([])]
        ok_code = bool(0)                             
            
    return mean_value,std_value,max_value,min_value,ok_code    
    
##########################################################################################
#useful founctions
##########################################################################################

def write_text_table_KROUPA(d,keys,fields,units,tab):
    #write to text table fields for given keys

    if len(units)>0:
        if not (len(fields)==len(units)):
            print '      *write_text_table_KROUPA - fields = ',fields,' - units = ',units,' - does not have the same length'
    print ' '
    
    #header
    line  = str('').ljust(tab[0])+'|'
    for i in range(0,len(fields)):
        string  = str(fields[i])
        if len(string)>tab[0]:
            string  = string[0:tab[0]-4]+' ...'
        line  = line + ' ' + string.ljust(tab[0])+'|'
    print line

    #units
    line  = str('').ljust(tab[0])+'|'
    for i in range(0,len(fields)):
        if (len(fields)==len(units)):
            string  = '['+str(units[i])+']'
        else:
            string  = '[-]'
        if len(string)>tab[0]:
            string  = string[0:tab[0]-4]+' ...'
        line  = line + ' ' + string.ljust(tab[0])+'|'
    print line
    
    #data
    for key in keys:
        line  = str(key).ljust(tab[0])+'|'
        for i in range(0,len(fields)):
            if check_existence_field_KROUPA(d[key],fields[i]):
                data_key,ok_code  = extract_data_field_KROUPA(d[key],fields[i])
                if len(units)>0 and len(data_key)>0:
                    if  ( (isinstance(units[0], int)) or (isinstance(units[0], float)) ) and ( (isinstance(data_key[0], int)) or (isinstance(data_key[0], float)) ):
                        number  = data_key[0]/units[i]
                        string  = str(number)
                    else: 
                        string  = str(data_key)
                else:
                    string  = str(data_key) 
            else:                
                string  = str('')
            if len(string)>tab[0]:
                string  = string[0:tab[0]-4]+' ...'
            line  = line + ' ' + string.ljust(tab[0])+'|'
        print line

    print ' '

def print_xls_file_KROUPA(d,id_file,max_len_string):    
    
    import xlwt
    
    known_fields  = []
    for key in d.keys():
        known_fields,ok_code_merge        = merge_lists_KROUPA(known_fields,d[key].keys())

    book = xlwt.Workbook(encoding="utf-8")
    sheet1 = book.add_sheet("data")

    #experiments
    column  = 0
    for field in known_fields:
        column  = column  + 1      
        string  = str(field)
        if len(string)>max_len_string:
            string  = string[0,max_len_string]
        sheet1.write(0, column, string)

    #data
    raw     =  0
    for key in d.keys():
        raw = raw + 1
        string  = str(key)
        if len(string)>max_len_string:
            string  = string[0,max_len_string]        
        sheet1.write(raw, 0, str(key))
        column  =   0
        for field in known_fields:
            column  = column  + 1      
            if check_existence_field_KROUPA(d[key],[field]):
                if not isinstance(d[key][field], dict):
                    string  = str(d[key][field])
                else:
                    string  = str('dict()')
                if len(string)>max_len_string:
                    string  = string[0:max_len_string]
                sheet1.write(raw, column, string)                

    book.save(id_file)
    
def synchronize_data_with_photographs_KROUPA(data,list_PHOTO,start_point_DATA,step_DATA):
    #Prepare data according the information given by step of photographs and start of the photographs and end of the photographs
    #input
    #data             = numpy array - data for resampling
    #list_PHOTO       = sorted list of raw photographs
    #start_PHOTO      = name of first photo to use
    #end_PHOTO        = name of last photo photo to use
    #start_point_DATA = value of starting point in data
    #step_DATA        = step of data
    #output
    #data_sync        = Prepared data in times when the photographs were taken
    #list_PHOTO_sync  = list of synchronized photographs
    #ok_code          = boolean true or false if everything goes well or not

    import numpy as np

    ok_code = bool(1)
    
    print '            CAUTION      *synchronize_data_with_photographs_KROUPA - In photo directory must be only pictures supposed to be used (and no subdirectories)!     CAUTION'
    print '                   list_PHOTO  = ',list_PHOTO
    print '                   start_POINT = ',start_point_DATA
    print '                   step_DATA   = ',step_DATA
    
    if ok_code:
    
        if start_point_DATA[0]<data[0]:
            print '      *synchronize_data_with_photographs_KROUPA - start point of data begins before data starts!'
    
        if start_point_DATA[0]>data[-1]:
            print '      *synchronize_data_with_photographs_KROUPA - start point of data begins after data ends!'
            ok_code = bool(0)
    
    data_sync       = np.array([])
    list_PHOTO_sync = []    
    
    if ok_code:        
        
        for index in range(0,len(list_PHOTO)):
        
            what_to_append  = start_point_DATA[0] + step_DATA[0] * float(index) 
            
            if what_to_append>=min(data) and what_to_append<=max(data):
            
                data_sync = np.append(data_sync,what_to_append)
            
                list_PHOTO_sync = list_PHOTO_sync + [list_PHOTO[index]]
                
            else:

                if ok_code:
                    if what_to_append>data[-1]:
                          print '      *synchronize_data_with_photographs_KROUPA - stop of appending data and photographs, because there are no data for next photographs!'
                          break
                
    return data_sync, list_PHOTO_sync, ok_code                       
    
    

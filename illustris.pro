; Illustris Simulation: Public Data Release.
; (util.pro): Various helper functions, and load all functions.

function partTypeNum, PT
  ; Mapping between common names and numeric particle types.
  compile_opt idl2, hidden, strictarr, strictarrsubs

  ; check if input is already numeric
  on_ioerror, check_names
  test = fix(PT)
  return, test
  
  check_names:
  partType = strlowcase(PT)

  if total(partType eq ['gas','cells']) then $
    return, 0
  if total(partType eq ['dm','darkmatter']) then $
    return, 1
  if total(partType eq ['tracer','tracers','tracermc','trmc']) then $
    return, 3
  if total(partType eq ['star','stars','stellar']) then $
    return, 4 ; only those with GFM_StellarFormationTime>0
  if total(partType eq ['wind']) then $
    return, 4 ; only those with GFM_StellarFormationTime<0
  if total(partType eq ['bh','bhs','blackhole','blackholes']) then $
    return, 5

  message,'Unknown particle type name.'
end

function str, tString
  compile_opt idl2, hidden, strictarr, strictarrsubs
  return, strcompress(string(tString),/remove_all)
end

function hdf5_all_attrs, file_obj, group_name
  ; return a hash() of all attributes within the group
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  g = h5g_open(file_obj, group_name)
  result = hash()
  
  for i=0,h5a_get_num_attrs(g)-1 do begin
    a_id = h5a_open_idx(g,i)
      result[h5a_get_name(a_id)] = h5a_read(a_id)
    h5a_close, a_id
  endfor
  
  h5g_close, g  
  return, result
end

function hdf5_dset_names, file_obj, group_name
  ; return a list() of all dataset names within the group
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  fields = list()
  for i=0,h5g_get_nmembers(file_obj,group_name)-1 do $
    fields.add, h5g_get_member_name(file_obj,group_name,i)
  return, fields
end

function hdf5_dset_properties, file_obj, group_name, shapes=shapes, types=types
  ; return a list() of data names in this group, optionally returning hashes via argument:
  ;  shapes: a {k:v} pair for each dataset in the group, where
  ;    k is the dataset name
  ;    v is a vector containing the dataspace size (one integer per dimension)
  ;  types: a {k:v} pair for each dataset in the group
  ;    k is the dataset name
  ;    v is the data type (?)
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  shapes = hash()
  types  = hash()  
  fields = list()
  
  g = h5g_open(file_obj, group_name)
  
  for i=0,h5g_get_nmembers(file_obj,group_name)-1 do begin
    dataset_name = h5g_get_member_name(file_obj,group_name,i)
    fields.add, dataset_name
    
    ; shape and type
    d = h5d_open(g,dataset_name)
      s = h5d_get_space(d)
        dataset_size = h5s_get_simple_extent_dims(s)
      h5s_close, s
      
      t = h5d_get_type(d)
        dataset_type = h5t_idltype(t)
      h5t_close, t
    h5d_close, d
    
    shapes[dataset_name] = dataset_size
    types[dataset_name]  = dataset_type
  endfor
  
  h5g_close, g
  
  return, fields
end

function hdf5_read_dataset_slice, file_obj, group_dataset_path, start, length
  ; read a hyperslab selection from a dataset
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if size(start,/dim) eq 0 then start = [start]
  if size(length,/dim) eq 0 then length = [length]
  
  d = h5d_open(file_obj, group_dataset_path)
    s = h5d_get_space(d)
      ; get dataspace size, prepare hyperslab indices
      dataset_size = h5s_get_simple_extent_dims(s)
      num_dims = n_elements(dataset_size)
      
      if size(start,/dim) ne num_dims or size(length,/dim) ne num_dims then $
        message,'Error: Input start and length vectors have incorrect dimensionality.'
      
      ; select hyperslab and create a memory space to hold the result (otherwise it is full size+sparse)
      h5s_select_hyperslab, s, start, length, /reset
      m = h5s_create_simple(length)
      
      ; read the data in the selected hyperslab
      data = h5d_read(d, file_space=s, memory_space=m)
      data = reform(data) ; remove any degenerate leading dimensions of size one
        
    h5s_close, s
  h5d_close, d
  
  return, data
end

; load other files
@illustris_groupcat
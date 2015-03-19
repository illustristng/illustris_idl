; Illustris Simulation: Public Data Release.
; groupcat.pro: File I/O related to the FoF and Subfind group catalogs.

function gcPath, basePath, snapNum, chunkNum=cn
  ; Return absolute path to a group catalog HDF5 file (modify as needed).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(cn) eq 0 then cn = 0
  
  gcPath = basePath + '/groups_' + string(snapNum,format='(I03)') + '/'
  filePath = gcPath + 'groups_' + string(snapNum,format='(I03)')
  filePath += '.' + str(cn) + '.hdf5'
  
  return, filePath  
end

function loadObjects, basePath, snapNum, gName, nName, fields=fields
  ; Load either halo or subhalo information from the group catalog.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  result = hash()
  
  ; load header from first chunk
  f = h5f_open( gcPath(basePath,snapNum) )
  header = hdf5_all_attrs(f, "Header")
  result['count'] = header['N'+nName+'_Total']
  
  if ~result['count'] then begin
    print, 'warning: zero groups, empty return (snap='+str(snapNum)+')'
    return, result
  endif
  
  ; if fields specified, verify existence
  field_names = hdf5_dset_properties(f, gName, shapes=shapes, types=types)
  
  if n_elements(fields) ne 0 then begin
    foreach field,fields do $
      if ~field_names.count(field) then $
        message,'Group catalog does not have requested field ['+field+']'
  endif
  
  ; if fields not specified, load everything
  if n_elements(fields) eq 0 then fields = field_names
  
  ; loop over all requested fields
  foreach field,fields do begin
    ; replace local length with global
    shape = shapes[field]
    shape[-1] = result['count']
    
    ; allocate within return hash
    result[field] = make_array(shape, type=types[field])
    
  endforeach
  
  h5f_close, f
  
  ; loop over chunks
  wOffset = 0
  
  for i=0,header['NumFiles']-1 do begin
    ; open chunk, load header and properties of datasets
    f = h5f_open( gcPath(basePath,snapNum,chunkNum=i) )
    
    header = hdf5_all_attrs(f, "Header")
    field_names = hdf5_dset_properties(f, gName, shapes=shapes, types=types)
    
    if ~header['N'+nName+'_ThisFile'] then continue ; empty file chunk
    
    ; loop over each requested field
    foreach field,fields do begin
      ; read data local to the current file
      shape = shapes[field]
      data = hdf5_read_dataset_slice(f, gName+"/"+field, lonarr(n_elements(shape)), shape)
      
      if n_elements(shape) eq 1 then $
        result[field,wOffset:wOffset+shape[-1]-1] = data
      if n_elements(shape) eq 2 then $
        result[field,*,wOffset:wOffset+shape[-1]-1] = data
          
    endforeach
    
    wOffset += shape[-1]
    
    h5f_close, f    
  endfor
  
  ; only a single field? then return the array instead of a single item hash
  if n_elements(fields) eq 1 then return, result[fields[0]]
  
  return, result
end

function loadSubhalos, basePath, snapNum, fields=fields
  ; Load all subhalo information from the entire group catalog for one snapshot
  ; (optionally restrict to a subset given by fields).
  return, loadObjects(basePath,snapNum,"Subhalo","subgroups",fields=fields)
end

function loadHalos, basePath, snapNum, fields=fields
  ; Load all halo information from the entire group catalog for one snapshot
  ; (optionally restrict to a subset given by fields).
  return, loadObjects(basePath,snapNum,"Group","groups",fields=fields)
end

function loadHeader, basePath, snapNum, chunkNum=cn
  ; Load the group catalog header.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  if n_elements(cn) eq 0 then cn = 0
  
  f = h5f_open( gcPath(basePath,snapNum,chunkNum=cn) )
  header = hdf5_all_attrs(f, "Header")
  h5f_close, f
  
  return, header
end
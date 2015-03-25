; Illustris Simulation: Public Data Release.
; illustris_groupcat.pro: File I/O related to the FoF and Subfind group catalogs.

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
  
  ; load header from first chunk
  f = h5f_open( gcPath(basePath,snapNum) )
  header = hdf5_all_attrs(f, "Header")
  
  result = hash()
  result['count'] = header['N'+nName+'_Total']
  
  if ~result['count'] then begin
    print, 'warning: zero groups, empty return (snap='+str(snapNum)+')'
    return, result
  endif
  
  ; if fields not specified, load everything
  field_names = hdf5_dset_properties(f, gName, shapes=shapes, types=types)
  
  if n_elements(fields) eq 0 then fields = field_names
  
  ; loop over all requested fields
  foreach field,fields do begin
    ; verify existence
    if ~field_names.count(field) then $
      message,'Group catalog does not have requested field ['+field+']'
      
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
      length = shapes[field]
      start  = lonarr(n_elements(length))
      data   = hdf5_read_dataset_slice(f, gName+"/"+field, start, length)
      
      ; save
      if n_elements(length) eq 1 then $
        result[field,wOffset:wOffset+length[-1]-1] = data
      if n_elements(length) eq 2 then $
        result[field,*,wOffset:wOffset+length[-1]-1] = data
    endforeach
    
    wOffset += length[-1]
    
    h5f_close, f    
  endfor
  
  ; only a single field? then return the array instead of a single item hash
  if n_elements(fields) eq 1 then return, result[fields[0]]
  
  return, result
end

function loadSubhalos, basePath, snapNum, fields=fields
  ; Load all subhalo information from the entire group catalog for one snapshot
  ; (optionally restrict to a subset given by fields).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  return, loadObjects(basePath,snapNum,"Subhalo","subgroups",fields=fields)
end

function loadHalos, basePath, snapNum, fields=fields
  ; Load all halo information from the entire group catalog for one snapshot
  ; (optionally restrict to a subset given by fields).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
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

function loadGroupcat, basePath, snapNum
  ; Load complete group catalog all at once.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  r = hash()
  r['subhalos'] = loadSubhalos(basePath, snapNum)
  r['halos']    = loadHalos(basePath, snapNum)
  r['header']   = loadHeader(basePath, snapNum)
  return, r
end

function loadGroupcatSingle, basePath, snapNum, haloID=hID, subhaloID=shID
  ; Return complete group catalog information for one halo or subhalo.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if (n_elements(hID) eq 0 and n_elements(shID) eq 0) or (n_elements(hID) and n_elements(shID)) then $
    message,'Error: Must specify either haloID or subhaloID (and not both).'

  if n_elements(hID)  gt 0 then gName = "Subhalo"
  if n_elements(shID) gt 0 then gName = "Group"
  if n_elements(hID)  gt 0 then searchID = hID
  if n_elements(shID) gt 0 then searchID = shID
  
  ; load groupcat offsets, calculate target file and offset
  f = h5f_open( gcPath(basePath,snapNum) )
    header = hdf5_all_attrs(f, "Header")
  h5f_close, f
  
  offsets = searchID - header['FileOffsets_'+gName]
  fileNum = max( where(offsets ge 0) )
  groupOffset = offsets[fileNum]
  
  ; load halo/subhalo fields into a hash
  result = hash()
  
  f = h5f_open( gcPath(basePath,snapNum,chunkNum=fileNum) )
  
    field_names = hdf5_dset_properties(f, gName, shapes=shapes)
    foreach field,field_names do begin
      ; parameters to read entire field of this chunk
      length = shapes[field]
      start  = lonarr(n_elements(length))
      
      ; modify to single element (could be multidimensional)
      start[-1]  = groupOffset
      length[-1] = 1
      
      result[field] = hdf5_read_dataset_slice(f, gName+"/"+field, start, length)
    endforeach
    
  h5f_close, f

  return, result
end

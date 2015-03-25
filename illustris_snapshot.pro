; Illustris Simulation: Public Data Release.
; illustris_snapshot.pro: File I/O related to the snapshot files.

function snapPath, basePath, snapNum, chunkNum=cn
  ; Return absolute path to a snapshot HDF5 file (modify as needed).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if n_elements(cn) eq 0 then cn = 0
  
  snapPath = basePath + '/snapdir_' + string(snapNum,format='(I03)') + '/'
  filePath = snapPath + 'snap_' + string(snapNum,format='(I03)')
  filePath += '.' + str(cn) + '.hdf5'
  
  return, filePath  
end

function getNumPart, header
  ; Calculate number of particles of all types given a snapshot header.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  nTypes = 6
  nPart  = ulon64arr(nTypes)
  
  for i=0,nTypes-1 do begin
    low_word  = ulong64(header['NumPart_Total',i])
    high_word = ulong64(header['NumPart_Total_HighWord',i])
    
    nPart[i] = ishft(high_word, 32) OR (low_word)
  endfor
  
  return, nPart
end

function loadSnapSubset, basePath, snapNum, partType, fields=fields, subset=subset
  ; Load a subset of fields for all particles/cells of a given partType.
  ; If offset and length specified, load only that subset of the partType.
  compile_opt idl2, hidden, strictarr, strictarrsubs  

  ptNum = partTypeNum(partType)
  gName = "PartType" + str(ptNum)
  
  ; load header from first chunk
  f = h5f_open( snapPath(basePath,snapNum) )
  header = hdf5_all_attrs(f, "Header")
  
  result = hash()
  nPart = getNumPart(header)
  
  ; decide global read size, starting file chunk, and starting file chunk offset
  if n_elements(subset) gt 0 then begin
    offsetsThisType = subset['offsetType',ptNum] - subset['snapOffsets',*,ptNum]
    
    fileNum = max(where( offsetsThisType ge 0 ))
    fileOff = offsetsThisType[fileNum]
    numToRead = subset['lenType',ptNum]
  endif else begin
    fileNum = 0
    fileOff = 0
    numToRead = nPart[ptNum]
  endelse
  
  result['count'] = numToRead  
  if ~numToRead then begin
    print, 'warning: no particles of requested type, empty return.'
    return, result
  endif
  
  ; find a chunk with this particle type
  i = 1
  while ~(hdf5_dset_names(f, "/")).count(gName) do begin
    h5f_close, f
    f = h5f_open( snapPath(basePath,snapNum,chunkNum=i) )
    i += 1
  endwhile
  
  ; if fields not specified, load everything
  field_names = hdf5_dset_properties(f, gName, shapes=shapes, types=types)
  
  if n_elements(fields) eq 0 then fields = field_names  

  ; loop over all requested fields
  foreach field,fields do begin
    print,field
    ; verify existence
    if ~field_names.count(field) then $
      message,'Particle type ['+str(ptNum)+'] does not have field ['+field+']'
      
    ; replace local length with global
    shape = shapes[field]
    shape[-1] = numToRead
    
    ; allocate within return hash
    result[field] = make_array(shape, type=types[field])
  endforeach
  
  h5f_close, f
  
  ; loop over chunks
  wOffset = 0
  origNumToRead = numToRead
  
  while numToRead gt 0 do begin
    print,fileNum,numToRead
    f = h5f_open( snapPath(basePath,snapNum,chunkNum=fileNum) )
    
    header = hdf5_all_attrs(f, "Header")
    field_names = hdf5_dset_properties(f, gName, shapes=shapes, types=types)
    
    ; no particles of requested type in this file chunk?
    if ~(hdf5_dset_names(f, "/")).count(gName) then begin
      if n_elements(subset) gt 0 then $
        message,'Read error: subset read should be contiguous.'
      fileNum += 1
      continue
    endif
    
    ; set local read length for this file chunk, truncate to be within the local size
    numTypeLocal = header['NumPart_ThisFile',ptNum]
    
    numToReadLocal = numToRead
    
    if fileOff + numToReadLocal gt numTypeLocal then $
      numToReadLocal = numTypeLocal - fileOff
      
    print,'['+string(fileNum,format='(I3)')+'] off='+str(fileOff)+' read ['+str(numToReadLocal)+$
          '] of ['+str(numTypeLocal)+'] remaining = '+str(numToRead-numToReadLocal)
    
    ; loop over each requested field for this particle type
    foreach field,fields do begin
      ; read data local to the current file
      length = shapes[field]
      start  = lonarr(n_elements(length))
      
      start[-1]  = fileOff
      length[-1] = numToReadLocal
      
      data = hdf5_read_dataset_slice(f, gName+"/"+field, start, length)
      
      ; save
      if n_elements(length) eq 1 then $
        result[field,wOffset:wOffset+length[-1]-1] = data
      if n_elements(length) eq 2 then $
        result[field,*,wOffset:wOffset+length[-1]-1] = data
    endforeach
    
    wOffset   += numToReadLocal
    numToRead -= numToReadLocal
    fileNum   += 1
    fileOff    = 0 ; start at beginning of all file chunks other than the first
    
    h5f_close, f
  endwhile
  
  if origNumToRead ne wOffset then $
    message,'Read ['+str(wOffset)+'] particles, but was expecting ['+str(origNumToRead)+']'

  ; only a single field? then return the array instead of a single item hash
  if n_elements(fields) eq 1 then return, result[fields[0]]
  
  return, result  
end

function getSnapOffsets, basePath, snapNum, id, type
  ; Compute offsets within snapshot for a particular group/subgroup.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  r = hash()
  
  ; load groupcat chunk offsets from header of first file
  f = h5f_open( gcPath(basePath,snapNum) )
    header = hdf5_all_attrs(f, "Header")
  h5f_close, f
  
  r['snapOffsets'] = header['FileOffsets_Snap']
  
  ; calculate target groups file chunk which contains this id
  groupFileOffsets = long(id) - header['FileOffsets_'+type]
  fileNum = max( where(groupFileOffsets ge 0) )
  groupOffset = groupFileOffsets[fileNum]
  
  ; load the length (by type) of this group/subgroup from the group catalog
  f = h5f_open( gcPath(basePath,snapNum,chunkNum=fileNum) )
    ; parameters to read entire field of this chunk
    field_names = hdf5_dset_properties(f, type, shapes=shapes)
    length = shapes[type+'LenType']
    start  = lonarr(n_elements(length))
    
    ; modify to single element (is multidimensional, [6,1])
    start[-1]  = groupOffset
    length[-1] = 1
    
    r['lenType'] = hdf5_read_dataset_slice(f, type+"/"+type+"LenType", start, length)
    
    ; load the offset (by type) of this group/subgroup within the snapshot
    r['offsetType'] = hdf5_read_dataset_slice(f, "Offsets/"+type+"_SnapByType", start, length)
  h5f_close, f
  
  return, r
end

function loadSubhalo, basePath, snapNum, id, partType, fields=fields
  ; load subhalo length, compute offset, call loadSnapSubset
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  subset = getSnapOffsets(basePath,snapNum,id,"Subhalo")
  return, loadSnapSubset(basePath,snapNum,partType,fields=fields,subset=subset)
end

function loadHalo, basePath, snapNum, id, partType, fields=fields
  ; load halo length, compute offset, call loadSnapSubset
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  subset = getSnapOffsets(basePath,snapNum,id,"Group")
  return, loadSnapSubset(basePath,snapNum,partType,fields=fields,subset=subset)
end

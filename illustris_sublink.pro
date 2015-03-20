; Illustris Simulation: Public Data Release.
; snapshot.pro: File I/O related to the Sublink merger tree files.

function sublinkPath, basePath, chunkNum=cn
  ; Return absolute path to a SubLink HDF5 file (modify as needed).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if n_elements(cn) eq 0 then cn = 0
  filePath = basePath + '/trees/SubLink/tree_extended.' + str(cn) + '.hdf5'  
  return, filePath  
end

function sublinkOffsets, basePath, snapNum, id
  ; Handle offset loading for a SubLink merger tree cutout.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  r = hash()
  
  ; load groupcat chunk offsets from header of first file
  f = h5f_open( gcPath(basePath,snapNum) )
    header = hdf5_all_attrs(f, "Header")
  h5f_close, f
  
  ; calculate target groups file chunk which contains this id
  groupFileOffsets = long(id) - header['FileOffsets_Subhalo']
  fileNum = max( where(groupFileOffsets ge 0) )
  groupOffset = groupFileOffsets[fileNum]
  
  ; load the length (by type) of this group/subgroup from the group catalog
  prefix = "Offsets/Subhalo_Sublink"
  fields = ['RowNum','LastProgenitorID','SubhaloID']
  
  f = h5f_open( gcPath(basePath,snapNum,chunkNum=fileNum) )
    ; read each single value (1d)
    foreach field,fields do $
      r[field] = hdf5_read_dataset_slice(f, prefix+field, groupOffset, 1)
  h5f_close, f
  
  return, r
end

function loadSublinkTree, basePath, snapNum, id, fields=fields, onlyMPB=onlyMPB
  ; Load portion of Sublink tree, for a given subhalo, in its existing flat format.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; the tree is all subhalos between SubhaloID and LastProgenitorID
  treeOff = sublinkOffsets(basePath,snapNum,id)
  
  rowStart = treeOff['RowNum']
  rowEnd   = treeOff['RowNum'] + (treeOff['LastProgenitorID'] - treeOff['SubhaloID'])
  nRows    = rowEnd - rowStart + 1

  ; find the tree file chunk containing this row
  f = h5f_open( gcPath(basePath,snapNum) )
    header = hdf5_all_attrs(f, "Header")
  h5f_close, f
  
  rowOffsets = rowStart - header['FileOffsets_SubLink',*]
  fileNum = max(where( rowOffsets ge 0 ))
  fileOff = rowOffsets[fileNum]
  
  ; load only main progenitor branch? in this case, get MainLeafProgenitorID now
  if keyword_set(onlyMPB) then begin
    print,fileNum,fileOff
    f = h5f_open( sublinkPath(basePath,chunkNum=fileNum) )
      MainLeafProgenitorID = hdf5_read_dataset_slice(f, 'MainLeafProgenitorID', fileOff, 1)
    h5f_close, f
    
    ; re-calculate nRows
    rowEnd = treeOff['RowNum'] + (MainLeafProgenitorID - treeOff['SubhaloID'])
    nRows  = rowEnd - rowStart + 1
  endif
  
  ; read
  result = hash()
  result['count'] = nRows
  
  f = h5f_open( sublinkPath(basePath,chunkNum=fileNum) )
    ; if no fields requested, return all fields
    field_names = hdf5_dset_properties(f, "/", shapes=shapes, types=types)
    
    if n_elements(fields) eq 0 then fields = field_names
    
    if fileOff + nRows gt (shapes['SubfindID'])[0] then $
      message,'Should not occur. Each tree is contained within a single file.'
      
    ; loop over each requested field
    foreach field,fields do begin
      ; verify existence
      if ~field_names.count(field) then $
        message,'SubLink tree does not have field ['+field+']'
        
      ; read
      length = shapes[field]
      start  = lonarr(n_elements(length))
      
      start[-1]  = fileOff
      length[-1] = nRows
      
      data = hdf5_read_dataset_slice(f, field, start, length)
      result[field] = data
    endforeach
  h5f_close, f
  
  ; only a single field? then return the array instead of a single item hash
  if n_elements(fields) eq 1 then return, result[fields[0]]
  
  return, result
end

function maxPastMass, tree, index, ptNum
  ; Get maximum past mass (of the given partType number) along the main branch of a subhalo
  ; specified by index within this tree.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  branchSize = tree['MainLeafProgenitorID',index] - tree['SubhaloID',index] + 1
  masses = tree['SubhaloMassType', ptNum, index : index + branchSize - 1]
  return, max(masses)
end

function numMergers, tree, minMassRatio=mMR, massPartType=mPT, index=ind
  ; Calculate the number of mergers in this sub-tree (optionally above some mass ratio threshold).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  ; set defaults
  if n_elements(mMR) eq 0 then mMR = 1e-10
  if n_elements(mPT) eq 0 then mPT = 'stars'
  if n_elements(index) eq 0 then index = 0
  
  ; verify the input sub-tree has the required fields
  reqFields = ['SubhaloID','NextProgenitorID','MainLeafProgenitorID','FirstProgenitorID','SubhaloMassType']
  
  foreach reqField,reqFields do $
    if ~tree.haskey(reqField) then $
      message, 'Error: Input tree needs to have loaded fields: '+strjoin(reqFields,",")

  numMergers = 0
  invMassRatio = 1.0 / mmR
  massPtNum = partTypeNum(mPT)
  
  ; walk back main progenitor branch
  rootID = tree['SubhaloID',index]
  fpID   = tree['FirstProgenitorID',index]
  
  while fpID ne -1 do begin
    fpIndex = index + (fpID - rootID)
    fpMass  = maxPastMass(tree, fpIndex, massPtNum)
    
    ; explore breadth
    npID = tree['NextProgenitorID',fpIndex]
  
    while npID ne -1 do begin
      npIndex = index + (npID - rootID)
      npMass  = maxPastMass(tree, npIndex, massPtNum)
      
      ; count if both masses are non-zero, and ratio exceeds threshold
      if fpMass gt 0.0 and npMass gt 0.0 then begin
        ratio = npMass / fpMass
        
        if ratio ge mMR and ratio le invMassRatio then numMergers += 1
      endif
      
      npID = tree['NextProgenitorID',npIndex]
    endwhile
    
    fpID = tree['FirstProgenitorID',fpIndex]
  endwhile
  
  return, numMergers  
end
      
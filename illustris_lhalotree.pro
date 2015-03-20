; Illustris Simulation: Public Data Release.
; illustris_lhalotree.pro: File I/O related to the LHaloTree merger tree files.

function lhalotreePath, basePath, chunkNum=cn
  ; Return absolute path to a LHaloTree HDF5 file (modify as needed).
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if n_elements(cn) eq 0 then cn = 0
  filePath = basePath + '/trees/treedata/trees_sf1_135.' + str(cn) + '.hdf5'  
  return, filePath  
end

function singleNodeFlat, conn, index, data_in, data_out, count, onlyMPB, gdp
  ; Recursive helper function: Add a single tree node.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  forward_function recProgenitorFlat
  
  if n_elements(data_in) gt 1 then begin
    ; in-memory walk (data_in is the already read data array)
    if size(data_out,/n_dim) eq 1 then $
      data_out[count] = data_in[index]
    if size(data_out,/n_dim) eq 2 then $
      data_out[*,count] = data_in[*,index]
  endif else begin
    ; on disk walk, read element we need (data_in is a hdf5 file_obj)
    start = gdp['start']
    start[-1] = index
    
    if size(data_out,/n_dim) eq 1 then $
      data_out[count] = hdf5_read_dataset_slice(data_in, gdp['dataset'], start, gdp['length'])
    if size(data_out,/n_dim) eq 2 then $
      data_out[*,count] = hdf5_read_dataset_slice(data_in, gdp['dataset'], start, gdp['length'])
  endelse
  
  count += 1
  
  count = recProgenitorFlat(conn,index,data_in,data_out,count,onlyMPB,gdp)
  return,count  
end

function recProgenitorFlat, conn, start_index, data_in, data_out, count, onlyMPB, gdp
  ; Recursive helper function: Flatten out the unordered LHaloTree, one data field at a time.
  compile_opt idl2, hidden, strictarr, strictarrsubs

  firstProg = conn['FirstProgenitor',start_index]
  
  if firstProg lt 0 then return, count
  
  ; depth-ordered traversal (down mpb)
  count = singleNodeFlat(conn,firstProg,data_in,data_out,count,onlyMPB,gdp)
  
  ; explore breadth
  if ~onlyMPB then begin
    nextProg = conn['NextProgenitor',firstProg]
    
    while nextProg ge 0 do begin
      count = singleNodeFlat(conn,nextProg,data_in,data_out,count,onlyMPB,gdp)
      nextProg = conn['NextProgenitor',nextProg]
    endwhile
  
  endif
  
  firstProg = conn['FirstProgenitor',firstProg]
  
  return,count
end
  
function loadLHaloTree, basePath, snapNum, id, fields=fields, onlyMPB=onlyMPB
  ; Load portion of LHaloTree, for a given subhalo, re-arranging into a flat format.
  compile_opt idl2, hidden, strictarr, strictarrsubs
  
  if n_elements(onlyMPB) eq 0 then onlyMPB = 0
  
  ; get offsets
  prefix    = "Offsets/Subhalo_LHalo"
  offFields = ['TreeFile','TreeIndex','TreeNum']
  treeOff   = treeOffsets(basePath,snapNum,id,prefix,offFields)
  
  ; config
  gName = 'Tree' + str(treeOff['TreeNum']) ; group name containing this subhalo
  nRows = !NULL ; we do not know in advance the size of the tree
  
  ; open
  fTree = h5f_open( lhalotreePath(basePath,chunkNum=treeOff['TreeFile']) )
    ; if no fields requested, return all fields
    field_names = hdf5_dset_properties(fTree, gName, shapes=shapes, types=types)
  
    if n_elements(fields) eq 0 then fields = field_names
    
    ; verify existence, loop over each requested field
    foreach field,fields do $
      if ~field_names.count(field) then message,'SubLink tree does not have field ['+field+']'

    ; load connectivity for this entire TreeX group
    connFields = ['FirstProgenitor','NextProgenitor']
    conn = hash()
    
    foreach field,connFields do $
      conn[field] = hdf5_read_dataset_slice(fTree, gName+"/"+field, 0, shapes[field])

    ; determine sub-tree size with dummy walk
    dummy = intarr(n_elements(conn['FirstProgenitor']))
    nRows = singleNodeFlat(conn, treeOff['TreeIndex'], dummy, dummy, 0, onlyMPB, 0)
    
    result = hash()
    result['count'] = nRows
    print, nRows
    
    ; walk through connectivity, one data field at a time
    foreach field,fields do begin
      print,field
      ; allocate the data array in the sub-tree
      length = shapes[field]
      start  = lonarr(n_elements(length))
      
      length[-1] = nRows
      data = make_array(length, type=types[field])
      
      ; load field for entire tree? doing so is much faster than randomly accessing the disk 
      ; during walk, assuming that the sub-tree is a large fraction of the full tree, and that 
      ; the sub-tree is large in the absolute sense. the decision is heuristic, and can be 
      ; modified (if you have the tree on a fast SSD, could disable the full load).
      
      if nRows lt 1000 then begin
        ; do not load, walk with single disk reads
        length[-1] = 1
        gdp = hash('dataset',gName+"/"+field, 'start',start, 'length',length)
        
        count = singleNodeFlat(conn, treeOff['TreeIndex'], fTree, data, 0, onlyMPB, gdp)
      endif else begin
        ; pre-load all, walk in-memory (use unmodified length)
        full_data = hdf5_read_dataset_slice(fTree, gName+"/"+field, start, shapes[field])

        count = singleNodeFlat(conn, treeOff['TreeIndex'], full_data, data, 0, onlyMPB, 0)
      endelse

      ; save field
      result[field] = data      
    endforeach    
  h5f_close, fTree
  
  ; only a single field? then return the array instead of a single item hash
  if n_elements(fields) eq 1 then return, result[fields[0]]
  
  return, result
end
    
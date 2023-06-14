; for THESAN only

FUNCTION cartPath, basePath, cartNum, chunkNum=cn
  ; Return absolute path to a cartesian HDF5 file (modify as needed).
  
  if n_elements(cn) eq 0 then cn = 0
  
  filePath_list = [basePath + '/cartesian_' + STRING(cartNum, FORMAT='(I03)') + '/cartesian_' + STRING(cartNum, FORMAT='(I03)') + '.' + str(cn) + '.hdf5']

  FOREACH filePath,filePath_list DO BEGIN
    IF FILE_TEST(filePath) THEN RETURN, filePath
  ENDFOREACH
  
  MESSAGE, 'No cartesian file found!'
  
END

FUNCTION getNumPixel, header
  ; Calculate number of pixels (per dimension) given a cartesian header.
  RETURN, long64(header['NumPixels'])
END

FUNCTION loadCartSubset, basePath, cartNum, fields=fields, bbox=bbox, sq=sq
  ; Load a subset of fields in the cartesian grids.
  ; If bbox is specified, load only that subset of data. bbox should have the 
  ; form [[start_i, start_j, start_k], [end_i, end_j, end_k]], where i,j,k are 
  ; the indices for x,y,z dimensions. Notice the last index is *inclusive*.
  ; If sq is True, return a numpy array instead of a dict if len(fields)==1.

  if n_elements(sq) eq 0 then sq = 1 ; by default true

  result = {}

  ; make sure fields is not a single element
  IF SIZE(fields, /TNAME) EQ 'STRING' THEN fields = [fields]

  ; load header from first chunk
  f = h5f_open( cartPath(basePath,cartNum,chunkNum=0) )
  header = hdf5_all_attrs(f, "Header")
  
  result = hash()
  nPix = getNumPixel(header)

  ; decide global read size, starting file chunk, and starting file chunk offset
  IF N_ELEMENTS(bbox) NE 0 THEN BEGIN
    load_all = 0
    start_i = bbox[0, 0]
    start_j = bbox[1, 0]
    start_k = bbox[2, 0]
    end_i = bbox[0, 1]
    end_j = bbox[1, 1]
    end_k = bbox[2, 1]
    IF start_i LT 0 THEN MESSAGE, 'start_i must be >= 0'
    IF start_j LT 0 THEN MESSAGE, 'start_j must be >= 0'
    IF start_k LT 0 THEN MESSAGE, 'start_k must be >= 0'
    IF end_i GE nPix THEN MESSAGE, 'end_i must be < nPix'
    IF end_j GE nPix THEN MESSAGE, 'end_j must be < nPix'
    IF end_k GE nPix THEN MESSAGE, 'end_k must be < nPix'
  END ELSE BEGIN
    load_all = 1
    bbox = [[0, 0, 0], [nPix-1, nPix-1, nPix-1]]
  END

  numToRead = long64(bbox[0,1]-bbox[0,0]+1) * long64(bbox[1,1]-bbox[1,0]+1) * long64(bbox[2,1]-bbox[2,0]+1)

  IF numToRead EQ 0 THEN RETURN, result

  ; if fields not specified, load everything
  field_names = hdf5_dset_properties(f, '/', shapes=shapes, types=types)
  
  IF n_elements(fields) eq 0 THEN BEGIN
      fields = field_names
  END

  ; loop over all requested fields
  FOREACH field,fields DO BEGIN
    ; verify existence
    if ~field_names.count(field) then $
      message,'Cartesian output does not have field ['+field+']'
      
    ; replace local length with global
    shape = shapes[field]
    shape[-1] = numToRead
    
    ; allocate within return hash
    result[field] = make_array(shape, type=types[field])

    print, field, shape

  ENDFOREACH
  
  h5f_close, f

  ; loop over chunks
  wOffset = 0l
  fileOffset = 0l
  origNumToRead = numToRead
  fileNum = 0l

  WHILE (numToRead GT 0) DO BEGIN

    f = h5f_open(cartPath(basePath, cartNum, chunkNum=fileNum))
    ; set local read length for this file chunk, truncate to be within the local size
    field_names = hdf5_dset_properties(f, '/', shapes=shapes, types=types)
    numPixelsLocal = shapes[field_names[0]]

    IF load_all THEN BEGIN
      pixToReadLocal = REPLICATE(1, numPixelsLocal)
      numToReadLocal = numPixelsLocal
    END ELSE BEGIN
      local_pixels_index = L64INDGEN(numPixelsLocal, start= fileOffset)
      local_pixels_i = local_pixels_index / (nPix^2)
      local_pixels_j = (local_pixels_index-local_pixels_i*nPix^2) / nPix
      local_pixels_k = local_pixels_index-local_pixels_i*nPix^2-local_pixels_j*nPix

      pixToReadLocal = WHERE((local_pixels_i GE bbox[0,0]) AND (local_pixels_i LE bbox[0,1]) AND $
                             (local_pixels_j GE bbox[1,0]) AND (local_pixels_j LE bbox[1,1]) AND $
                             (local_pixels_k GE bbox[2,0]) AND (local_pixels_k LE bbox[2,1]), numToReadLocal)
    END

    ; loop over each requested field for this particle type
    FOREACH field,fields DO BEGIN
      length = shapes[field]
      start  = lonarr(n_elements(length))
      data = hdf5_read_dataset_slice(f, "/"+field,start,length)
      result[field,wOffset:wOffset+numToReadLocal-1] = data[pixToReadLocal]
    ENDFOREACH

    wOffset   += numToReadLocal
    numToRead -= numToReadLocal

    fileOffset += numPixelsLocal
    fileNum += 1

    print, fileNum, wOffset, numToRead, numToReadLocal 
    
    h5f_close, f
  ENDWHILE

  ; verify we read the correct number
  IF origNumToRead NE wOffset THEN MESSAGE, 'Read [' + STRING(wOffset) + '] particles, but was expecting [' + STRING(origNumToRead) + ']'

  ; only a single field? then return the array instead of a single item dict
  IF sq AND (N_ELEMENTS(fields) EQ 1) THEN RETURN, result[0]
  RETURN, result
END

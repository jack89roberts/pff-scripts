function [ FONTData ] = loadExtractedFONTData( fileName )
%loadExtractedFONTData Loads an extracted FONT data file and returns the
%data entered in to the FONTData struct.
        FONTData = load(fileName);
        FONTData = FONTData.FONTData;
end


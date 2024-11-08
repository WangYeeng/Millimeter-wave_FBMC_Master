function CodingMatrix = GetPrecodingMatrix(obj,StartIndex)
     % Returns a precoding matrix which allows QAM transmission in FBMC-OQAM.
     % The second argument represents the start index  of the 
     % Precoding matrix and can be set to 1 or 2. 
     
     if rem(log2(obj.Nr.MCSymbols),1)
        error('Number of time-symbols must be a power of 2');
     end              
     PrecodingMatrix = fwht(eye(obj.Nr.MCSymbols),obj.Nr.MCSymbols,'sequency')*sqrt(obj.Nr.MCSymbols);
     CodingMatrix_Basis(:,:,1) = PrecodingMatrix(:,1:2:end);
     CodingMatrix_Basis(:,:,2) = PrecodingMatrix(:,2:2:end);

     IndexSpreading = repmat((1:obj.Nr.Subcarriers)',[1,obj.Nr.MCSymbols]);
     IndexSpreadingHalf = IndexSpreading(:,1:2:end);
     CodingMatrix = zeros(obj.Nr.Subcarriers*obj.Nr.MCSymbols,obj.Nr.Subcarriers*obj.Nr.MCSymbols/2);
     for i_CodeMapping = 1:max(max(IndexSpreading))
         CodingMatrix(IndexSpreading(:)==i_CodeMapping,IndexSpreadingHalf(:)==i_CodeMapping) = ...
         CodingMatrix_Basis(:,:,mod(i_CodeMapping+StartIndex,2)+1);        
     end    
             
end  
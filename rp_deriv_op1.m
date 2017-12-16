function [M1a,M1b] = rp_deriv_op1(data,deriv_flag)
%% set up 2-D derivative operator

% data is the input matrix

% allocate space for the sparse operator
M1a = spalloc(size(data,1)-1,size(data,1)-1,(2*(size(data,1)-1)));

M1b = spalloc(size(data,2)-1,size(data,2)-1,(2*(size(data,2)-1)));

totalcounty=1;
totalcountx=1;


switch deriv_flag
    case 1
        % first derivative operator
        
        for count=1:size(data,2)
            
            for counter=1:size(data,1)-1
                
                diaglocx=(count-1)*size(data,1)+counter;
                diaglocy=totalcounty;
                
                M1a(diaglocy,diaglocx)=-1;
                M1a(diaglocy,diaglocx+1)=1;
                
                
                totalcounty=totalcounty+1;
            end
            
        end
        
    case 2
        for count=1:size(data,2)
            
            for counter=1:size(data,1)-2
                
                diaglocx=(count-1)*size(data,1)+counter;
                diaglocy=totalcounty;
                
                M1a(diaglocy,diaglocx)=1;
                M1a(diaglocy,diaglocx+1)=-2;
                M1a(diaglocy,diaglocx+2)=1;
                
                
                totalcounty=totalcounty+1;
            end
            
        end
end


%%
totalcounty=1;
totalcountx=1;

switch deriv_flag
    case 1
        for count=1:size(data,2)-1
            
            for counter=1:size(data,1)
                
                diaglocx=(count-1)*size(data,1)+counter;
                diaglocy=totalcountx;
                
                M1b(diaglocy,diaglocx)=-1;
                M1b(diaglocy,diaglocx+size(data,1))=1;
                
                totalcountx=totalcountx+1;
            end
            
        end
        
    case 2
        for count=1:size(data,2)-2
            
            for counter=1:size(data,1)
                
                diaglocx=(count-1)*size(data,1)+counter;
                diaglocy=totalcountx;
                
                M1b(diaglocy,diaglocx)=1;
                M1b(diaglocy,diaglocx+size(data,1))=-2;
                M1b(diaglocy,diaglocx+(2*size(data,1)))=1;
                
                totalcountx=totalcountx+1;
            end
            
        end
end

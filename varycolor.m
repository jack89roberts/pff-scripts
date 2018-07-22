function ColorSet=varycolor(NumberOfPlots)
% VARYCOLOR Produces colors with maximum variation on plots with multiple
% lines.
%
%     VARYCOLOR(X) returns a matrix of dimension X by 3.  The matrix may be
%     used in conjunction with the plot command option 'color' to vary the
%     color of lines.  
%
%     Yellow and White colors were not used because of their poor
%     translation to presentations.
% 
%     Example Usage:
%         NumberOfPlots=50;
%
%         ColorSet=varycolor(NumberOfPlots);
% 
%         figure
%         hold on;
% 
%         for m=1:NumberOfPlots
%             plot(ones(20,1)*m,'Color',ColorSet(m,:))
%         end

%Created by Daniel Helmick 8/12/2008

narginchk(1,1)%correct number of input arguements??
nargoutchk(0, 1)%correct number of output arguements??

%Take care of the anomolies. Edited by Jack to give default colours for
%MatLab 2014b+
if NumberOfPlots<1
    ColorSet=[];

elseif NumberOfPlots<7
    ColorSet=zeros(NumberOfPlots,3);
    ColorSet(1,:)=[0, 0.4470, 0.7410];
   
    if NumberOfPlots>=2
        ColorSet(2,:)=[0.8500, 0.3250, 0.0980];
    end
    if NumberOfPlots>=3
        ColorSet(3,:)=[0.9290, 0.6940, 0.1250];
    end
    if NumberOfPlots>=4
        ColorSet(4,:)=[0.4940, 0.1840, 0.5560];
    end
    if NumberOfPlots>=5
        ColorSet(5,:)=[0.4660, 0.6740, 0.1880];
    end
    if NumberOfPlots>=6
        ColorSet(6,:)=[0.3010, 0.7450, 0.9330];
    end   
    
else %default and where this function has an actual advantage

    %we have 5 segments to distribute the plots
    EachSec=floor(NumberOfPlots/5); 
    
    %how many extra lines are there? 
    ExtraPlots=mod(NumberOfPlots,5); 
    
    %initialize our vector
    ColorSet=zeros(NumberOfPlots,3);
    
    %This is to deal with the extra plots that don't fit nicely into the
    %segments
    Adjust=zeros(1,5);
    for m=1:ExtraPlots
        Adjust(m)=1;
    end
    
    SecOne   =EachSec+Adjust(1);
    SecTwo   =EachSec+Adjust(2);
    SecThree =EachSec+Adjust(3);
    SecFour  =EachSec+Adjust(4);
    SecFive  =EachSec;

    for m=1:SecOne
        ColorSet(m,:)=[0 1 (m-1)/(SecOne-1)];
    end

    for m=1:SecTwo
        ColorSet(m+SecOne,:)=[0 (SecTwo-m)/(SecTwo) 1];
    end
    
    for m=1:SecThree
        ColorSet(m+SecOne+SecTwo,:)=[(m)/(SecThree) 0 1];
    end
    
    for m=1:SecFour
        ColorSet(m+SecOne+SecTwo+SecThree,:)=[1 0 (SecFour-m)/(SecFour)];
    end

    for m=1:SecFive
        ColorSet(m+SecOne+SecTwo+SecThree+SecFour,:)=[(SecFive-m)/(SecFive) 0 0];
    end
    
end
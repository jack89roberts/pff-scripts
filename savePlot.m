function [ status ] = savePlot( saveDir, saveName, saveFormats, overwrite, quitOnError )
%savePlot Saves current figure in EPS, PDF, PNG and FIG formats.
%
%   status = savePlot( saveDir, saveName, saveFormats, overwrite, quitOnError )
%
%   saveDir - directory to save figure in. If it doesn't exist it is
%   created. Directory name should not have trailing /.
%
%   saveName - file name for plot, to which the fig/png/eps extension will
%   be added.
%
%   overwrite - If 1 will overwrite pre-existing plots at same location,
%   if 0 adds time stamp to saveName if there is a conflict, if -1 doesn't
%   save any new files if something already exists with same name 
%   (default: 1).
%
%   saveFormats - [x,x,x,x], x=0 or 1. First element=.png., second=.fig,
%   third=.eps, fourth=.pdf. Default: [1,1,0,1] (save .fig, .png and .pdf)
%
%   quitOnError - if true quits if there's a problem, if false
%   just gives a warning (default: false).
%
%   status: 000 - no issues, 100 - directory issue, 010 - file
%   overwrite warning, 001 - save files issue. Plus combinations. 

    if (nargin < 5 || isempty(quitOnError))
        quitOnError = 0;
    end
    if (nargin < 4 || isempty(overwrite))
        overwrite = 1;
    end
    if (nargin < 3 || isempty(saveFormats))
        saveFormats = [1 1 0 1];
    end
    if (nargin<2 || isempty(saveName))
        % assume directory and name without extension given together as
        % first argument
        [saveDir,saveName] = fileparts(saveDir);
        if (isempty(saveDir) || isempty(saveName))
            error('savePlot:saveDir','Invalid input format for saveDir');
        end
    end
    if (nargin<1 || isempty(saveDir))
        error('savePngEpsFig:nargin','Incorrect number of input arguments');
    end
    
    status = [0,0,0];
    
    % if the save directory doesn't exist, create it
    if (~isdir(saveDir))
        try 
            mkdir(saveDir);
        catch
            status(1) = 1;
            if (quitOnError)
                error('savePngEpsFig:mkdir','Error making directory %s',saveDir);
            else
                warning('savePngEpsFig:mkdir','Error making directory %s',saveDir);
            end    
        end
    end
    
    % construct file path
    saveStr = [saveDir '/' saveName];
   
    % check if files already exist
    if (overwrite ~= 1)
        if ( (exist([saveStr '.eps'],'file') == 2) ||  (exist([saveStr '.png'],'file') == 2) ||  (exist([saveStr '.fig'],'file') == 2) )
            status(2) = 1;

    %         if(overwrite == 1) % just give a warning that some files will be overwritten
    %             warning('savePngEpsFig:overwrite','Files with name like %s will be overwritten',saveStr);

            if(overwrite == -1) % return from function as user said not to save any new files if they already exist
                warning('savePngEpsFig:overwrite','No new files saved there are already file(s) with name %s',saveStr);
                return;

            else % overwrite == 0 or some weird value, take safest option and save files with time stamp added to avoid overwriting anything     
                timeNow = datevec(now);
                yearNow = timeNow(1);
                monthNow = timeNow(2);
                dayNow = timeNow(3);
                hourNow = timeNow(4);
                minuteNow = timeNow(5);
                secondNow = round(timeNow(6));

                saveStr = sprintf('%s_%04d%02d%02d_%02d%02d%02d',saveStr,yearNow,monthNow,dayNow,hourNow,minuteNow,secondNow);

                warning('savePngEpsFig:overwrite','Time stamp added to new files with name like %s to avoid overwriting old files',saveStr);
            end

        end
    end
    
    try
        if (saveFormats(1))
            print([saveStr '.png'],'-dpng');
        end
        if (saveFormats(2))
            savefig([saveStr '.fig']);
        end
        if (saveFormats(3))
            print([saveStr '.eps'],'-depsc');
        end
        if (saveFormats(4))
            h = gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);       
            print([saveStr '.pdf'], '-dpdf');
        end
    catch
        status(3) = 1;
        if (quitOnError)
            error('savePngEpsFig:save','Error saving file(s) with names like %s',saveStr);
        else
            warning('savePngEpsFig:save','Error saving file(s) with names like %s',saveStr);
        end    
    end

end


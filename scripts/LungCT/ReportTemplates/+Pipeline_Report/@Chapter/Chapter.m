classdef Chapter < mlreportgen.report.Chapter 

    properties 
        DATE
    end 

    methods 
        function obj = Chapter(varargin) 
            obj = obj@mlreportgen.report.Chapter(varargin{:}); 
        end 
    end 

    methods (Hidden) 
        function templatePath = getDefaultTemplatePath(~, rpt) 
            path = abc.Chapter.getClassFolder(); 
            templatePath = ... 
                mlreportgen.report.ReportForm.getFormTemplatePath(... 
                path, rpt.Type); 
        end 

    end 

    methods (Static) 
        function path = getClassFolder() 
            [path] = fileparts(mfilename('fullpath')); 
        end 

        function createTemplate(templatePath, type) 
            path = abc.Chapter.getClassFolder(); 
            mlreportgen.report.ReportForm.createFormTemplate(... 
                templatePath, type, path); 
        end 

        function customizeReporter(toClasspath) 
            mlreportgen.report.ReportForm.customizeClass(... 
                toClasspath, "abc.Chapter"); 
        end 

    end  
end
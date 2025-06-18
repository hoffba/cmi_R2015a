classdef MimServer < handle
    
    properties (Access = private)
        WebSocketServer
        Mim
        ModelController
        Reporting
    end
    
    methods (Static)
        function start()
            MimServer.getServer().startServer();
        end
        
        function stop()
            MimServer.getServer().stopServer();
        end
        
        function clear()
            MimServer.getServer().clearModels();
        end
        
        function localUpdateModel(modelName, hash, data)
            MimServer.getServer().updateLocalModelValue(modelName, hash, data);            
        end
        
        function triggerLocalUpdate(modelName)
            MimServer.getServer().triggerUpdateLocalModelValue(modelName);            
        end
        
        
        
        function mimServer = getServer()
            persistent mimServerSingleton
            if isempty(mimServerSingleton) || ~isvalid(mimServerSingleton)
                mimServerSingleton = MimServer();
            end
            mimServer = mimServerSingleton;            
        end
    end
    
    methods (Access = private)
        function obj = MimServer()
            framework_def = PTKFrameworkAppDef;
            obj.Reporting = MimReporting([], [], 'mimserver.log');
            obj.Mim = MimMain(framework_def, obj.Reporting);
            obj.ModelController = MimModelController(obj.Mim);
            obj.WebSocketServer = MimWebSocketServer(30000, obj.ModelController);
        end
        
        function startServer(obj)
            if ~obj.WebSocketServer.Status
                obj.WebSocketServer.start();
            end
        end
        
        function stopServer(obj)
            obj.WebSocketServer.stop();
        end

        function clearModels(obj)
            obj.ModelController.clear();
            obj.WebSocketServer.clearModels();
        end
        
        function updateLocalModelValue(obj, modelName, hash, value)
            obj.WebSocketServer.updateLocalModelValue(modelName, hash, value);
        end
        
        function triggerUpdateLocalModelValue(obj, modelName)
            [value, hash] = obj.ModelController.getValue(modelName);
            obj.WebSocketServer.updateLocalModelValue(modelName, hash, value);
        end
        
    end
end

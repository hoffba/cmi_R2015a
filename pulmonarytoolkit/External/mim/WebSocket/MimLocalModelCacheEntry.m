classdef MimLocalModelCacheEntry < MimModelCacheEntry

    properties
        CachedValue
        Controller
        ModelId
    end
    
    events
        ValueChanged
    end
    
    properties (Access = private)
        Listeners
    end
    
    methods
        function obj = MimLocalModelCacheEntry(currentHash, remoteHash, controller, modelId)
            obj = obj@MimModelCacheEntry(currentHash, remoteHash);
            obj.Controller = controller;
            obj.ModelId = modelId;
            obj.Listeners = [];
        end
        
        function updateAndNotify(obj, currentHash, remoteHash, value)
            % Updates the hashes and value based on values from a remote source.
            % This notifies local listeners
            obj.updateHashes(currentHash, remoteHash);
            obj.CachedValue = value;
            
            % Notifies the local model framework
            obj.Controller.setValue(obj.ModelId, value);
        end
        
        function modifyCurrentHashAndValue(obj, newHash, newValue)
            % Updates the current hash and value based on values from a local source.
            % This does not update the RemoteHash since the source is
            % local. Local listeners are not notified since it is assumed
            % that they are already aware of local changes.
            obj.CurrentHash = newHash;
            obj.CachedValue = newValue;
        end
        
        function cache = getCache(obj)
        	cache = obj;
        end
        
        function value = getCurrentValue(obj)
            value = obj.CachedValue;
        end
    end
end


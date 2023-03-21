function opts = init_pipeline_opts

opts = struct( ...
              'dcm_path',       ''      ,... % Location of catalog file
              'save_path',      ''      ,... % Location to save pipeline results (for GL, must be on Turbo)
              'par_size',       5       ,... % parallel processes for local analysis
              'cluster',        'batch' ,... % Where to process the data: 'GL', 'batch', or 'debug'
              'cluster_type',   'auto'  ,... % GL cluster type to use: 'largemem', 'standard'
              'mem',            24      ,... % GL memory needed per process
              'nnodes',         1       ,... % GL number of nodes (only for standard)
              'username',       ''      ,... % GL your uniquename for email notifications
              ...
              'unreg',          true    ,... % MODULE unregistered stats for Exp/Ins
              'airway',         true    ,... % MODULE airways results from YACTA
              'scatnetAT',      true    ,... % MODULE scatternet air trapping
              'vessel',         true    ,... % MODULE blood vessel analysis
              'reg',            true    ,... % MODULE register Ins to Exp
              'scatnetEmph',    true    ,... % MODULE scatternet Emph map
              'prm',            true    ,... % MODULE PRM analysis
              'tprm',           true    ,... % MODULE tPRM analysis
              'dBlood',         true    ,... % MODULE generation of dBlood map
              ...
              'quickreg',       false   ,... % OPTION to skip last step in registration for speed
              'orient_check',   true    ,... % OPTION to auto-correct image orientation
              'jac',            true    ,... % OPTION to generate Jacobian from reg
              'jacmat',         false   ,... % OPTION to generate Jacobian matrix from reg
              'def',            false   ,... % OPTION to generate deformation fields from reg
              'Tvoi',           true     ... % OPTION to transform Ins VOI
             );

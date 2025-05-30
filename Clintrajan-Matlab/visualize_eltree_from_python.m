function visualize_eltree_from_python(tree_elpi, X, finalCTAtable, feature_name, output_path)
    % Create Python script for visualization

    % Determine relevant file names
    savepath = fileparts(output_path);
    fn_script = fullfile(savepath, 'temp_visualize_eltree.py');
    fn_input = fullfile(savepath, 'temp_viz_input.mat');
    fn_features = fullfile(savepath, 'temp_feature_name.txt');

    python_path = update_python_paths();
    python_script = sprintf([...
        'import numpy as np\n',...
        'import matplotlib.pyplot as plt\n',...
        'import networkx as nx\n',...
        'from scipy.io import savemat, loadmat\n',...
        'import os\n\n',...
        'def visualize_eltree_with_data(tree, X, color_values, feature_name, output_file):\n',...
        '    # Extract tree components\n',...
        '    nodep = tree["NodePositions"]\n',...
        '    edges = tree["Edges"]\n',...
        '    \n',...
        '    # Create graph\n',...
        '    g = nx.Graph()\n',...
        '    g.add_edges_from(edges)\n',...
        '    pos = nx.kamada_kawai_layout(g, scale=2)\n',...
        '    \n',...
        '    # Convert positions to array\n',...
        '    idx = np.array([pos[j] for j in range(len(pos))])\n',...
        '    \n',...
        '    # Partition data\n',...
        '    partition = np.zeros(len(X), dtype=int)\n',...
        '    dists = np.zeros(len(X))\n',...
        '    \n',...
        '    # Assign each point to nearest node\n',...
        '    for i in range(len(X)):\n',...
        '        dists_to_nodes = np.sum((nodep - X[i])**2, axis=1)\n',...
        '        partition[i] = np.argmin(dists_to_nodes)\n',...
        '        dists[i] = np.min(dists_to_nodes)\n',...
        '    \n',...
        '    # Project points onto graph\n',...
        '    projval = np.zeros(len(X))\n',...
        '    edgeid = np.zeros(len(X), dtype=int)\n',...
        '    X_proj = np.zeros_like(X)\n',...
        '    \n',...
        '    for i in range(len(X)):\n',...
        '        min_dist = float("inf")\n',...
        '        best_edge = 0\n',...
        '        best_t = 0\n',...
        '        best_proj = None\n',...
        '        \n',...
        '        for j in range(len(edges)):\n',...
        '            node1 = edges[j][0]\n',...
        '            node2 = edges[j][1]\n',...
        '            \n',...
        '            p1 = nodep[node1]\n',...
        '            p2 = nodep[node2]\n',...
        '            \n',...
        '            # Project point onto line segment\n',...
        '            v = p2 - p1\n',...
        '            u = X[i] - p1\n',...
        '            \n',...
        '            # Calculate projection parameter\n',...
        '            v_dot_v = np.dot(v, v)\n',...
        '            if v_dot_v > 0:  # Avoid division by zero\n',...
        '                t = max(0, min(1, np.dot(u, v) / v_dot_v))\n',...
        '            else:\n',...
        '                t = 0\n',...
        '            \n',...
        '            # Calculate projected point\n',...
        '            proj = p1 + t * v\n',...
        '            \n',...
        '            # Calculate distance\n',...
        '            d = np.sum((X[i] - proj)**2)\n',...
        '            \n',...
        '            if d < min_dist:\n',...
        '                min_dist = d\n',...
        '                best_edge = j\n',...
        '                best_t = t\n',...
        '                best_proj = proj\n',...
        '        \n',...
        '        projval[i] = best_t\n',...
        '        edgeid[i] = best_edge\n',...
        '        X_proj[i] = best_proj\n',...
        '    \n',...
        '    # Calculate distance to projection and add scatter\n',...
        '    dist2proj = np.sum(np.square(X - X_proj), axis=1)\n',...
        '    shift = np.percentile(dist2proj, 20)\n',...
        '    dist2proj = dist2proj - shift\n',...
        '    dist2proj[dist2proj < 0] = 0\n',...
        '    \n',...
        '    # Generate scattered points around tree\n',...
        '    scatter_parameter = 0.025\n',...
        '    x = np.zeros(len(X))\n',...
        '    y = np.zeros(len(X))\n',...
        '    \n',...
        '    for i in range(len(X)):\n',...
        '        # Distance from edge\n',...
        '        r = 0\n',...
        '        if dist2proj[i] > 0:\n',...
        '            r = np.sqrt(dist2proj[i]) * scatter_parameter\n',...
        '        \n',...
        '        # Get node coordinates for this edge\n',...
        '        if edgeid[i] < len(edges):\n',...
        '            edge = edges[edgeid[i]]\n',...
        '            x_coos = [idx[edge[0]][0], idx[edge[1]][0]]\n',...
        '            y_coos = [idx[edge[0]][1], idx[edge[1]][1]]\n',...
        '            \n',...
        '            projected_on_edge = False\n',...
        '            \n',...
        '            if projval[i] < 0:\n',...
        '                # Project to 0%% of edge\n',...
        '                x_coo = x_coos[0]\n',...
        '                y_coo = y_coos[0]\n',...
        '            elif projval[i] > 1:\n',...
        '                # Project to 100%% of edge\n',...
        '                x_coo = x_coos[1]\n',...
        '                y_coo = y_coos[1]\n',...
        '            else:\n',...
        '                # Project to appropriate %% of edge\n',...
        '                x_coo = x_coos[0] + (x_coos[1] - x_coos[0]) * projval[i]\n',...
        '                y_coo = y_coos[0] + (y_coos[1] - y_coos[0]) * projval[i]\n',...
        '                projected_on_edge = True\n',...
        '            \n',...
        '            # Calculate vector orthogonal to edge\n',...
        '            vex = x_coos[1] - x_coos[0]\n',...
        '            vey = y_coos[1] - y_coos[0]\n',...
        '            \n',...
        '            if not projected_on_edge:\n',...
        '                vex = np.random.random() - 0.5\n',...
        '                vey = np.random.random() - 0.5\n',...
        '            \n',...
        '            vn = np.sqrt(vex*vex + vey*vey)\n',...
        '            if vn > 0:\n',...
        '                vex = vex/vn\n',...
        '                vey = vey/vn\n',...
        '                rsgn = 1 if np.random.random() < 0.5 else -1\n',...
        '                x[i] = x_coo + vey * r * rsgn\n',...
        '                y[i] = y_coo - vex * r * rsgn\n',...
        '            else:\n',...
        '                x[i] = x_coo\n',...
        '                y[i] = y_coo\n',...
        '        else:\n',...
        '            x[i] = 0\n',...
        '            y[i] = 0\n',...
        '    \n',...
        '    # Create figure\n',...
        '    plt.figure(figsize=(12, 9))\n',...
        '    plt.style.use("ggplot")\n',...
        '    \n',...
        '    # Scatter points\n',...
        '    plt.scatter(x, y, c=color_values, cmap="jet", s=2, alpha=0.7)\n',...
        '    \n',...
        '    # Add colorbar\n',...
        '    plt.colorbar(label=feature_name)\n',...
        '    \n',...
        '    # Scatter nodes\n',...
        '    plt.scatter(idx[:,0], idx[:,1], s=10, c="black", alpha=0.8)\n',...
        '    \n',...
        '    # Plot edges\n',...
        '    for j in range(len(edges)):\n',...
        '        x_coo = [idx[edges[j][0]][0], idx[edges[j][1]][0]]\n',...
        '        y_coo = [idx[edges[j][0]][1], idx[edges[j][1]][1]]\n',...
        '        plt.plot(x_coo, y_coo, c="k", linewidth=1.5, alpha=0.5)\n',...
        '    \n',...
        '    # Set title\n',...
        '    plt.title(f"Tree Visualization with {feature_name}", fontsize=14)\n',...
        '    \n',...
        '    # Save figure\n',...
        '    plt.savefig(output_file, dpi=300, bbox_inches="tight")\n',...
        '    plt.close()\n\n',...
        'def process_tree_data(input_file, output_file):\n',...
        '    # Load data\n',...
        '    data = loadmat(input_file)\n',...
        '    X = data["X"]\n',...
        '    tree = {\n',...
        '        "NodePositions": data["NodePositions"],\n',...
        '        "Edges": data["Edges"]\n',...
        '    }\n',...
        '    color_values = data["color_values"].flatten()\n',...
        '    \n',...
        '    # Get feature name from text file\n',...
        '    with open("%s", "r") as f:\n',...
        '        feature_name = f.read().strip()\n',...
        '    \n',...
        '    print(f"Feature name: {feature_name}")\n',...
        '    \n',...
        '    # Create visualization\n',...
        '    visualize_eltree_with_data(tree, X, color_values, feature_name, output_file)\n',...
        '\n',...
        'if __name__ == "__main__":\n',...
        '    import sys\n',...
        '    process_tree_data(sys.argv[1], sys.argv[2])\n'],strrep(fn_features,'\','\\'));
    
    % Save Python script
    fid = fopen(fn_script, 'w');
    fprintf(fid, '%s', python_script);
    fclose(fid);
    
    % Find feature values
    feature_idx = find(strcmp(finalCTAtable.Properties.VariableNames, feature_name));
    
    if isempty(feature_idx)
        warning('Feature %s not found in finalCTAtable', feature_name);
        color_values = ones(size(X, 1), 1);
    else
        color_values = table2array(finalCTAtable(:, feature_idx(1)));
    end
    
    % Extract NodePositions and Edges from tree_elpi
    NodePositions = tree_elpi.NodePositions;
    Edges = tree_elpi.Edges;
    
    % Save feature name as a text file
    fid = fopen(fn_features, 'w');
    fprintf(fid, '%s', feature_name);
    fclose(fid);
    
    % Save variables for Python
    save(fn_input, 'X', 'NodePositions', 'Edges', 'color_values');
    
    % Prepare command with full paths
    command = sprintf('"%s" "%s" "%s" "%s"', ...
        python_path, ...
        fn_script, ...
        fn_input, ...
        output_path);
    
    % Execute Python script
    [status, commandOutput] = system(command);

    % Check for errors
    if status ~= 0
        fprintf('Error output from Python:\n%s\n', commandOutput);
        error('Python script failed to execute properly');
    else
        fprintf('%s\n', commandOutput);  % Display Python output
    end
    
    % Clean up temporary files
    delete(fn_input);
    delete(fn_script);
    delete(fn_features);
    
    fprintf('Tree visualization saved to %s\n', output_path);
end

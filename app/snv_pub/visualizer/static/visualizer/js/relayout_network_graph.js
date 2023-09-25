document
    .getElementById('relayout-btn')
    .addEventListener('click', async () => {
        const selectorPairs = [
            ['.layout-parameter[name="marker.size"]', '#traces-to-update-node-size input:checked[type="checkbox"]'],
            ['.layout-parameter[name="line.width"]', '#traces-to-update-line-width input:checked[type="checkbox"]']
        ]; 
        
        for (const [parameterSelector, traceSelector] of selectorPairs) {
            relayoutPlotly('network-graph', parameterSelector, traceSelector);
        }
        
    });

const relayoutPlotly = (graphId, parameterSelector, traceSelector='') => {
    console.log('Start relayoutPlotly()');
    const graph = document.getElementById(graphId);
    console.log(graph);
    if (graph.data === undefined) {
        console.log('Do nothing');
        ;
    } else {
        const parameters = document.querySelectorAll(parameterSelector);

        console.log(parameters);
        let data_update = {};
        let layout_update = {};
        for (const parameter of parameters) {
            if (parameter.value !== '') {
                console.log(parameter.className);
                console.log(typeof parameter.className);
                const parameterClasses = parameter.className.split(' ');
                console.log(parameterClasses);

                let parameterValue = parameter.value;
                if (parameter.type === 'number') {
                    parameterValue = Number(parameter.value);
                }

                if (parameterClasses.includes('data_update')) {
                    data_update[parameter.name] = parameterValue;
                } else if (parameterClasses.includes('layout_update')) {
                    layout_update[parameter.name] = parameterValue;
                } else {
                    ;
                }
            }
        }
        
        console.log(`traceSelector: ${traceSelector}`);
        const traces = document.querySelectorAll(traceSelector);
        console.log(traces);
        let traceIdxList = [];
        for (const trace of traces) {
            traceIdxList.push(trace.value)
        }
        
        console.log(data_update);
        console.log(layout_update);
        /* When traceIdxList is empty, all traces will be updated. */
        Plotly.update(graphId, data_update, layout_update, traceIdxList);
    }
}

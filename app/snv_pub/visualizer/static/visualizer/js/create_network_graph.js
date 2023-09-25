const getCookie = (name) => {
    if (document.cookie && document.cookie !== '') {
        for (const cookie of document.cookie.split(';')) {
            const [key, value] = cookie.trim().split('=')
            if (key === name) {
                return decodeURIComponent(value)
            }
        }
    }
}

const csrfToken = getCookie('csrftoken')


const removeAllChildren = (parentElement) => {
    return new Promise((resolve, reject) => {
        if (parentElement === undefined) {
            reject(`parantElement is undefined`);
        } else {
            while (true) {
                if (parentElement.firstChild) {
                    parentElement.removeChild(parentElement.firstChild);
                } else  {
                    resolve(parentElement);
                    break;
                }
            }
        }
    });
}


const displayMessage = (messageUlSelector, listMessageHTMLVsLevel) => {
    const messageUl = document.querySelector(messageUlSelector);
    messageUl.classList.add('messages');
    messageUl.classList.remove('empty');

    removeAllChildren(messageUl).then((messageUlCleaned) => {
        for (const [messageHTML, level] of listMessageHTMLVsLevel) {
            const messageLiHTML = `<li class="message ${level}">${messageHTML}</li>`
            messageUlCleaned.insertAdjacentHTML('beforeend', messageLiHTML);
        }
    })
}

const displayFormFieldValidation = (fieldErrors) => {
    for (const fieldError of fieldErrors) {
        console.log(fieldError);
        console.log(`#${fieldError.fieldId} + ul.errorlist`);
        const messageContainer = document.querySelector(`#${fieldError.fieldId} + ul.errorlist`);
        if (messageContainer) {
            console.log(messageContainer);
            removeAllChildren(messageContainer).then((messageContainerCleand) => {
                messageContainerCleand.classList.add(fieldError.level);
                messageContainerCleand.classList.remove('empty')

                for (const message of fieldError.messages) {
                    const messageHTML = `<li class="message ${fieldError.level}">${message.message}</li>`;
                    messageContainerCleand.insertAdjacentHTML('beforeend', messageHTML);
                }
            });            
        }

        // const messageContainer = document.getElementById(fieldError.messageContainerId);
        // if (messageContainer) {
        //     messageContainer.classList.add(fieldError.level);
        //     messageContainer.classList.remove('empty')
        //     const messageHTML = `<p class="message ${fieldError.level}">${fieldError.message}</p>`;
        //     messageContainer.insertAdjacentHTML('beforeend', messageHTML);
        // }
    }    
}

const removeAllFormFieldValidation = (selectorOfMessageContainer) => {
    const messageContainers = document.querySelectorAll(selectorOfMessageContainer);
    console.log(messageContainers);
    for (const messageContainer of messageContainers) {
        removeAllChildren(messageContainer);
    }
}


plotlyDatarevision = 0;
const createNetworkGraph = (formID, graphId) => {
    return new Promise((resolve, reject) => {
        console.time('Get response')
        const url = document.getElementById(formID).action;
        const inputAllFiles = document.querySelectorAll('input[type="file"]');
        const inputParameters = document.querySelectorAll('.network-parameter:not(input[type="file"]):not(input[name="adducts_for_suspect_compound"])');
        const adductsForSuspectCompound = document.querySelectorAll('input:checked[name="adducts_for_suspect_compound"]');
        const submittedButton = document.activeElement;

        let body = new FormData();

        inputAllFiles.forEach((inputFile) => {
            console.log(inputFile);
            console.log(inputFile.files[0]);
            for (const file of inputFile.files) {
                body.append(inputFile.name, file, file.name);
            }
        })

        body.append('create_or_update', submittedButton.value);

        // console.log(`Node indexes to remove: ${Array.from(setOfNodeIdxToRemove)}`);
        setOfNodeIdxToRemove.forEach((nodeIdx) => {
            body.append('l_total_input_idx_to_remove', nodeIdx);
        })        

        for (const parameter of inputParameters) {
            if (parameter.type === 'checkbox') {
                // console.log(`${parameter.checked}, ${typeof parameter.checked}`)
                body.append(parameter.name, parameter.checked);
            } else if (parameter.value !== '') {
                body.append(parameter.name, parameter.value);
            }
        }

        for (const adduct of adductsForSuspectCompound) {
            if (adduct.disabled === false) {
                body.append(adduct.name, adduct.value);
            }
        }

        // Delete all field validation errors.
        removeAllFormFieldValidation('.form-field-message-container:not(.empty)');

        // body.set('subgraph_depth', '12.5');
        // body.set('score_threshold', '-1');
        
        // body.set('quant_polar', 'aaa');
        // body.set('node_select_subgraph_mode', 'foo');
        // body.set('subgraph_num_core_nodes', 'aaa');
        // body.set('quant_value', 'aaa');
        // body.set('stat_value', 'aaa');
        // body.set('quant_subgraph_depth', 'aaa');

        // body.set('filter_select_category', 'aaa');
        // body.set('mz_tolerance', 'aaa');
        // body.set('ion_mode', 'aaa');
        // body.append('adducts_for_suspect_compound', 'foo');
        // body.append('adducts_for_suspect_compound', 'bar');
        // body.set('str_mass_defects', '0.1,aaa');
        console.log(body)
        fetch(url, {
            method: 'POST',
            body: body,
            headers: {
                'X-CSRFToken': csrfToken
            },
        }).then((response) => {
            console.timeEnd('Get response')
            console.log(response);
            console.log(response.headers.get('Content-Type'));
            if (!response.ok) {
                throw new Error(response.statusText);
            }
            return response;
        }).then((response) => {
            const content_type = response.headers.get('Content-Type');
            console.log(content_type);
            if (content_type.includes('application/json')) {
                return response.json();                    
            }
        }).then((responseJson) => {
            console.log(responseJson);
                
            if (!responseJson.formValidation.isValid) {
                const fieldErrors = responseJson.fieldErrors;
                displayFormFieldValidation(fieldErrors);
                throw new Error('Invalid input');
            }
            console.time('Create plot')
            const listLayerNameVsTrace = responseJson.listLayerNameVsTrace;
            // edgeTraceBetweenSampleAndRef is contained in listLayerNameVsTrace
            // const edgeTraceBetweenSampleAndRef = responseJson.edgeTraceBetweenSampleAndRef;
            const edgeTraceInterSample = responseJson.edgeTraceInterSample;
            const listMeshAnnotation = responseJson.listMeshAnnotation;
            const layout = responseJson.layout;

            let data = [];
            for (const [layerName, trace] of listLayerNameVsTrace) {
                data.push(trace);
            }
            // data.push(edgeTraceBetweenSampleAndRef, edgeTraceInterSample);
            data.push(edgeTraceInterSample);
            
            plotlyDatarevision += 1;
            layout.datarevision = plotlyDatarevision;
            
            Plotly.react(graphId, data, layout);
            console.timeEnd('Create plot')
            console.time('Customize plot')
            const networkPlot = document.getElementById(graphId);
            resizeObserverOfPlotly.observe(networkPlot);
            
            // Remove all event handlers from networkPlot
            networkPlot.removeAllListeners();
            networkPlot.on('plotly_click', (data) => {
                displayPointData('point-data-container', data);
                displayCompoundSVG('compound-structure-container', data);
            });

            // let legendClickCount = 0;
            // networkPlot.on('plotly_legendclick', (data) => {
            //     /*
            //     If this function returns false, default function for 'plotly_legendclick' event will not be executed
            //     and 'plotly_legenddoubleclick' event will never be fired.
            //     */
            //     ++legendClickCount;

            //     setTimeout(() => {
            //         if (legendClickCount === 1) {
            //             // Single-click
            //             customLegendClick(data, graphId);
            //         } else if (legendClickCount === 2) {
            //             // Double-click
            //             customLegendDoubleClick(data, graphId);
            //         }
                    
            //         // Initialize legendClickCount
            //         legendClickCount = 0;
            //     }, 350);

            //     return false;
            // });
            
            // Hide an element for Plotly's hover label.
            let defaultHoverLabel = document.querySelector(`#${graphId} #scene svg`)
            if (! defaultHoverLabel) {
                defaultHoverLabel = document.querySelector(`#${graphId} g.hoverlayer`);
            }
            defaultHoverLabel.style.display = 'none';

            console.log(networkPlot.data);

            networkPlot.on('plotly_hover', (data) => {
                displayCustomHover('my-hover-container', data);
            });

            networkPlot.on('plotly_unhover', (data) => {
                let myHover = document.getElementById('my-hover-container');
                myHover.style.visibility = 'hidden';
                // myHover.style.display = 'none';
            });
            
            console.timeEnd('Customize plot')
            displayMessage('.message-container .messages', []);

            resolve(data);

        }).catch((error) => {
            if (error.message) {
                displayMessage('.message-container .messages', [[error.message, 'error']]);
            }            
            console.error(error);
            reject(error);
        });
    });
}

const resizeObserverOfPlotly = new ResizeObserver((entries) => {
    for (const entry of entries) {
        if (entry.target.data !== undefined) {
            const currentWidth = entry.contentRect.width;
            const currentHeight = entry.contentRect.height;

            Plotly.relayout(entry.target.id, {width: currentWidth, height: currentHeight});
        }
    }
});

const getHoverTextAsHTML = (data) => {
    let text = `x: ${data.points[0].x}<br>y: ${data.points[0].y}<br>z: ${data.points[0].z}<br>`;

    if (data.points[0]['data']['mode'] === 'lines') {
        text = '';
    }

    if (data.points[0].text !== undefined) {
        text += `Attributes ------------------------------<br>${data.points[0].text}`;
    }
    if (data.points[0]['marker.color'] !== undefined) {
        text += `<br>Marker color: ${data.points[0]['marker.color']}`;
    }
    
    return text;
}

const getHoverTextAsString = (data) => {
    let text = `x: ${data.points[0].x}\ny: ${data.points[0].y}\nz: ${data.points[0].z}`;

    if (data.points[0].text !== undefined) {
        const hoverText = data.points[0].text.replaceAll('<br>', '\n');
        text += `\nAttributes: ------------------------------\n${hoverText}`;
    }
    if (data.points[0]['marker.color'] !== undefined) {
        text += `\nMarker color: ${data.points[0]['marker.color']}`;
    }
    
    return text;
}

/**
 * 
 * @param {String} hoverText - Hover text. (e.g. 'key1: value1<br>key2: value2<br>key3: value3')
 * @param {String} delimiter - Delimiter between key and value.
 * @return {Array.<Array>} [[String(key1), String(value1)], [String(key2), String(value2)], [String(key2), String(value2)]]
 */
const convertHoverTextToArray = (hoverText, delimiter) => {
    let hoverData = [];
    const lines = hoverText.split(/\n|(<br>)/);

    for (const line of lines) {
        const keyAndValueString = line.split(delimiter);

        if (keyAndValueString.length === 1) {
            ;
        } else if (keyAndValueString.length === 2) {
            let value = keyAndValueString[1].trim();
            if (value.match(/(nan)|(NaN)|(None)/)) {
                value = ''
            }
            hoverData.push([keyAndValueString[0], value]);
        } else {
            hoverData.push([keyAndValueString[0], keyAndValueString.slice(1).join(delimiter).trim()]);
        }
    }

    return hoverData;
}

const getGlobalAccessionFromHoverText = (hoverText, delimiter) => {
    const lines = hoverText.split(/\n|(<br>)/);
    let globalAccession = '';

    for (const line of lines) {
        const keyAndValueString = line.split(delimiter);
        
        if (keyAndValueString[0] === 'Global accession') {
            if (keyAndValueString.length === 1) {
                globalAccession = '';
            } else if (keyAndValueString.length === 2) {
                globalAccession = keyAndValueString[1].trim();
                if (globalAccession.match(/(nan)|(NaN)|(None)/)) {
                    globalAccession = ''
                }
            } else {
                globalAccession = keyAndValueString.slice(1).join(delimiter).trim()
            }
            break;
        }        
    }
    // console.log(globalAccession);
    return globalAccession;
}

const getGlobalAccessionArr = (traces) => {
    const arrHoverData = [];
    for (const trace of traces) {
        // console.log(trace);
        let accessions = [];
        if (trace.type === 'scatter3d') {
            if (trace.mode === 'markers') {
                accessions = trace.text.map((x) => getGlobalAccessionFromHoverText(x, ':'))
            }
        }
        arrHoverData.push(accessions);
    }
    // console.log(arrHoverData);
    return arrHoverData;
}

const getNodeColorArr = (traces) => {
    const arrNodeColor = [];
    for (const trace of traces) {
        // console.log(trace);
        let colors = [];
        if (trace.type === 'scatter3d') {
            if (trace.mode === 'markers') {
                colors = trace.marker.color;
            }
        }
        arrNodeColor.push(colors);
    }
    // console.log(arrNodeColor);
    return arrNodeColor;
}


const displayCustomHover = (parentId, data) => {
    const pageX = data.points[0].bbox.x0;
    const pageY = data.points[0].bbox.y0;
    let pointColor = '#000000';
    let textColor = '#ffffff';
    let svg_base64 = ''
    if (data.points[0].customdata) {
        if (data.points[0].customdata[0]) {
            pointColor = data.points[0].customdata[0];
        }
        if (data.points[0].customdata[1]) {
            textColor = data.points[0].customdata[1];
        }
        if (data.points[0].customdata[2]) {
            svg_base64 = data.points[0].customdata[2];
        }
    }

    const myHover = document.getElementById(parentId);
    removeAllChildren(myHover).then((myHoverCleaned) => {
        myHoverCleaned.style.top = pageY + 10; // 10px is padding of the <section id="main-network-graph-section"> element.
        myHoverCleaned.style.left = pageX + 10;
        myHoverCleaned.style['background-color'] = pointColor;
        myHoverCleaned.style.visibility = 'visible';
        // myHoverCleaned.style.display = '';

        const hoverText = getHoverTextAsHTML(data);
        let hoverTextHTML = '';
        if (svg_base64) {
            hoverTextHTML = `<img style="width: 200px; height: 200px; background-color: #ffffff;" src="${svg_base64}">\n`
        }
        hoverTextHTML += `<p class="custom-hover-text" style="background-color: transparent;">${hoverText}</p>`;
        hoverTextHTML = `<div class="custom-hover-content" style="background-color: ${pointColor}; color: ${textColor}; z-index: ${myHoverCleaned.style["z-index"] + 1};">${hoverTextHTML}</div>`
        myHoverCleaned.insertAdjacentHTML('afterbegin', hoverTextHTML);
    });
}

const displayPointData = (parentId, data) => {
    /* Display data of a clicked point */
    // console.log(data);
    const pointDataContainer = document.getElementById(parentId);
    const text = getHoverTextAsHTML(data);
    const hoverDataList = convertHoverTextToArray(text, ':');
    // console.log(hoverDataList);
    
    for (const hoverData of hoverDataList) {
        if (hoverData[0] === 'ID') {
            clickedNode.idx = hoverData[1];
            break;
        }
    }
    // console.log(clickedNode);

    let trHTML = '';
    for (const hoverData of hoverDataList) {
        trHTML += `<tr><th>${hoverData[0]}</th><td>${hoverData[1]}</td></tr>`;
    }
    const tableHTML = `<table class="table-vertical-headers"><tbody>${trHTML}</tbody></table>`;

    removeAllChildren(pointDataContainer).then((pointDataContainerCleaned) => {
        pointDataContainerCleaned.insertAdjacentHTML('beforeend', tableHTML);
    });
}

const displayCompoundSVG = (parentId, data) => {
    const compoundSVGContainer = document.getElementById(parentId);
    let svg_base64 = ''
    if (data.points[0].customdata) {
        if (data.points[0].customdata[0]) {
            pointColor = data.points[0].customdata[0];
        }
        if (data.points[0].customdata[1]) {
            textColor = data.points[0].customdata[1];
        }
        if (data.points[0].customdata[2]) {
            svg_base64 = data.points[0].customdata[2];
        }
    }

    removeAllChildren(compoundSVGContainer).then((compoundSVGContainerCleaned) => {
        if (svg_base64) {
            imageHTML = `<img src="${svg_base64}">\n`
            compoundSVGContainerCleaned.insertAdjacentHTML('beforeend', imageHTML);
        }
    });    
}

const getTraceIndexVsName = (data) => {
    let arrTraceNameVsIndex = [];

    if (data.length === 0) {
        ;
    } else {
        data.forEach((trace, idx) => {
            let traceName = trace.name;
            if (traceName === undefined) {
                traceName = '';
            }

            arrTraceNameVsIndex.push([idx, traceName]);
        });
    }
    return arrTraceNameVsIndex;
}

const createOptions = (selectId, options) => {
    const select = document.getElementById(selectId);

    if (select === null) {
        ;
    } else {
        /* Delete checkbox elements in parentElement. */
        const currentOptions = document.querySelectorAll(`#${selectId} option`);
        for (const currentOption of currentOptions) {
            currentOption.remove();
        }

        for (const [value, text] of options) {
            const uuid = generateUuid();
            optionHTML = `<option id="${uuid}" value="${value}">${text}</option>`
            select.insertAdjacentHTML('beforeend', optionHTML);
        }
    }
}

const createCheckboxes = (parentId, checkboxes) => {
    const parentElement = document.getElementById(parentId);

    if (parentElement === null) {
        ;
    } else {
        /* Delete checkbox elements in parentElement. */
        const currentChexboxes = document.querySelectorAll(`#${parentId} div.checkbox-label-container`);
        for (const currentChexbox of currentChexboxes) {
            currentChexbox.remove();
        }

        for (const [value, text] of checkboxes) {
            const uuid = generateUuid();
            const checkboxHTML = `<div class="flex-row checkbox-label-container">
                <input id="${uuid}" type="checkbox" value="${value}">
                <label for="${uuid}" class="flex-row">
                    ${text}
                    <span class="material-symbols-outlined blank-checkbox">check_box_outline_blank</span>
                    <span class="material-symbols-outlined selected-checkbox">select_check_box</span>
                </label>
            </div>`;
            parentElement.insertAdjacentHTML('beforeend', checkboxHTML);
            // console.log('Add checkbox');
        }
    }
}

const createCancelList = (parentId, cancelList, mySet = new Set()) => {
    const parentElement = document.getElementById(parentId);
    if (parentElement === null) {
        ;
    } else {
        /* Delete cancel-item-container elements in parentElement. */
        const currentCancelList = document.querySelectorAll(`#${parentId} div.cancel-item-container`);
        for (const currentCancelItem of currentCancelList) {
            currentCancelItem.remove();
        }

        for (const [value, text] of cancelList) {
            const uuid = generateUuid();
            const checkboxHTML = `<div id="container-${uuid}" class="flex-row cancel-item-container">
                <button id="cancel-btn-${uuid}" type="button" class="cancel-button" value="${value}" title="cancel">
                    <span class="material-icons" value="${value}">cancel</span>
                </button>
                <p id="${uuid}" value="${value}">${text}</p>
            </div>`;
            parentElement.insertAdjacentHTML('beforeend', checkboxHTML);

            document.getElementById(`cancel-btn-${uuid}`)
                .addEventListener('click', (e) => {
                    document.getElementById(`container-${uuid}`).remove();
                    if (mySet.size) {
                        mySet.delete(e.target.attributes.value.value);
                        // console.log(`mySet: ${Array.from(mySet)}`);
                    }
                });
        }
    }
}

const generateUuid = () => {
    let chars = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.split("");
    chars.forEach((char, idx) => {
        if (char === 'x') {
            chars[idx] = Math.floor(Math.random() * 16).toString(16);
        } else if (char === 'y') {
            chars[idx] = (Math.floor(Math.random() * 4) + 8).toString(16);
        } else {
            ;
        }
    });

    return chars.join('');
}




const customLegendClick = (data, graphId) => {
    /* Mesh3D layer will be hidden if its node and edge trace are hidden. */

    let layerToShow = new Set();

    const clickedTrace = data.data[data.curveNumber]
    let visibleOfClickedTrace = 'legendonly';
    if (clickedTrace.visible === 'legendonly') {
        visibleOfClickedTrace = true;
        layerToShow.add(clickedTrace.legendgroup)
    }

    data.data.forEach((trace, index) => {
        if (trace.type === 'mesh3d') {
            ;
        } else {
            if (trace.visible !== 'legendonly' && index !== data.curveNumber) {
                layerToShow.add(trace.legendgroup)
            }            
        }
    });

    let layerIndexToShow = [];
    let layerIndexToHide = [];
    data.data.forEach((trace, index) => {
        if (trace.type === 'mesh3d') {
            if (layerToShow.has(trace.legendgroup)) {
                layerIndexToShow.push(index);
            } else {
                layerIndexToHide.push(index);
            }
        }
    });

    Plotly.update(graphId, {visible: visibleOfClickedTrace}, {}, [data.curveNumber]);

    if (layerIndexToShow.length) {
        Plotly.update(graphId, {visible: true}, {}, layerIndexToShow);
    }
    if (layerIndexToHide.length) {
        Plotly.update(graphId, {visible: false}, {}, layerIndexToHide);
    }    
}

const customLegendDoubleClick = (data, graphId) => {
    /* Mesh3D layer containing clicked trace is will not be hidden. */
    const layerToShow = data.data[data.curveNumber].legendgroup;
    
    let traceIndexAll = [];
    let traceIndexesToShow = [];
    let traceIndexesToHide = [];

    data.data.forEach((trace, index) => {
        traceIndexAll.push(index);

        if (trace.legendgroup === layerToShow) {
            traceIndexesToShow.push(index);
        } else {
            traceIndexesToHide.push(index);
        }
    });

    for (const trace of data.data) {
        if (trace.visible === 'legendonly') {
            traceIndexesToShow = traceIndexAll;
            traceIndexesToHide = [];
            break;
        }
    }

    if (traceIndexesToShow.length) {
        Plotly.update(graphId, {visible: true}, {}, traceIndexesToShow);
    }
    if (traceIndexesToHide.length) {
        Plotly.update(graphId, {visible: 'legendonly'}, {}, traceIndexesToHide);
    }
    
}

document.addEventListener('DOMContentLoaded', () => {
    // 
    document.getElementById('spectral-network-form').addEventListener('submit', async (e) => {
        e.preventDefault();
        if (document.activeElement.id !== 'btn-to-create-network-graph'
            && document.activeElement.id !== 'btn-to-update-network-graph'
            && document.activeElement.id !== 'btn-to-execute-node-removal') {
            throw new Error();
        }

        document.getElementById('network-graph').style.display = 'none';
        document.getElementById('plotly-loading-img-container').style.display = '';
        createNetworkGraph('spectral-network-form', 'network-graph').then((figureData) => {
            document.getElementById('plotly-loading-img-container').style.display = 'none';
            const arrGlobalAccession = getGlobalAccessionArr(figureData);
            const arrNodeColor = getNodeColorArr(figureData);
            window.sessionStorage.setItem(['strArrGlobalAccession'], [JSON.stringify(arrGlobalAccession)]);
            window.sessionStorage.setItem(['strarrNodeColor'], [JSON.stringify(arrNodeColor)]);
            window.sessionStorage.setItem('strArrNodeColorbyTraceOriginal', null);
            window.sessionStorage.setItem('traceIdxChangedNodeColor', null);

            return getTraceIndexVsName(figureData);
        }).then((arrTraceIndexVsName) => {
            let nodeTraces = [];
            let edgeTraces = [];
            for (const [traceIdx, traceName] of arrTraceIndexVsName) {
                if (traceName.startsWith('Compounds of ')) {
                    nodeTraces.push([traceIdx, traceName])
                } else if (traceName.startsWith('Similarities between ')) {
                    edgeTraces.push([traceIdx, traceName]);
                } else {
                    ;
                }
            }

            document.querySelector('#graph-section .overlay').style.display = 'none';

            createCheckboxes('traces-to-update-node-size', nodeTraces);
            createCheckboxes('traces-to-update-line-width', edgeTraces);
            document.querySelector('#graph-layout-container .overlay').style.display = 'none';
            document.getElementById('network-graph').style.display = '';
        }).catch(() => {document.getElementById('plotly-loading-img-container').style.display = 'none';});
    });

    // Create a set of node indexes to remove.
    clickedNode = new Object();
    setOfNodeIdxToRemove = new Set();
    document.getElementById('btn-to-remove-node').addEventListener('click', () => {
        if(clickedNode.idx) {
            setOfNodeIdxToRemove.add(clickedNode.idx);
            
            let cancelList = [];
            setOfNodeIdxToRemove.forEach((nodeIdx) => {cancelList.push([nodeIdx, nodeIdx])});

            createCancelList('container-of-nodes-to-remove', cancelList, mySet=setOfNodeIdxToRemove);
        }
        // console.log(`setOfNodeIdxToRemove: ${Array.from(setOfNodeIdxToRemove)}`);
    })

    document.querySelector('#node-removal-section .cancel-all-button').addEventListener('click', () => {
        setOfNodeIdxToRemove.clear();
        createCancelList('container-of-nodes-to-remove', []);
    });

    // Change node color
    document.getElementById('btn-to-search-node').addEventListener('click', () => {
        const globalAccession = document.getElementById('input-node-to-search').value;
        const nodeColor = document.getElementById('node-color-to-search').value;
        let elemResult = document.getElementById('accession-search-result');
        let resultMessage = '';
        let fragGlobalAccessionMatch = false;

        if (globalAccession) {
            const strArrGlobalAccession = window.sessionStorage.getItem(['strArrGlobalAccession']);
            const arrGlobalAccession = JSON.parse(strArrGlobalAccession);
            const strarrNodeColor = window.sessionStorage.getItem(['strarrNodeColor']);
            let arrNodeColor = JSON.parse(strarrNodeColor);

            const strArrNodeColorbyTraceOriginal = window.sessionStorage.getItem('strArrNodeColorbyTraceOriginal');
            if (strArrNodeColorbyTraceOriginal !== 'null') {
                // console.log(`strArrNodeColorbyTraceOriginal: ${strArrNodeColorbyTraceOriginal}`);

                const arrNodeColorbyTraceOriginal = JSON.parse(strArrNodeColorbyTraceOriginal);
                // console.log(`arrNodeColorbyTraceOriginal: ${arrNodeColorbyTraceOriginal}`);

                let traceIdxChangedNodeColor = window.sessionStorage.getItem('traceIdxChangedNodeColor');
                
                traceIdxChangedNodeColor = Number(traceIdxChangedNodeColor);
                // console.log(`traceIdx: ${traceIdxChangedNodeColor}`);

                Plotly.restyle('network-graph', {'marker.color': [arrNodeColorbyTraceOriginal]}, [traceIdxChangedNodeColor]);
            }
            

            // console.log(arrGlobalAccession);
            let traceIdx = 0;
            for (const arrGlobalAccessionbyTrace of arrGlobalAccession) {
                if (!arrGlobalAccessionbyTrace) {
                    continue;
                }

                const accessionIdx = arrGlobalAccessionbyTrace.indexOf(globalAccession);

                if (accessionIdx !== -1) {
                    let arrNodeColorbyTrace = arrNodeColor[traceIdx];

                    window.sessionStorage.setItem('strArrNodeColorbyTraceOriginal', JSON.stringify(arrNodeColorbyTrace));
                    window.sessionStorage.setItem('traceIdxChangedNodeColor', traceIdx);

                    arrNodeColorbyTrace[accessionIdx] = nodeColor;
                    // console.log(arrNodeColorbyTrace);
                    // console.log(`traceIdx: ${traceIdx}`);
                    Plotly.restyle('network-graph', {'marker.color': [arrNodeColorbyTrace]}, [traceIdx]);
                    fragGlobalAccessionMatch = true;
                    resultMessage = 'Found!'
                    break;
                }
                traceIdx += 1;
            }
            
            if (!fragGlobalAccessionMatch) {
                resultMessage = `Not found: ${globalAccession}`
            }
            elemResult.innerText = resultMessage;
            elemResult.style.visibility = 'visible';
        }
        
    });
})


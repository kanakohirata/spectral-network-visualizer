document.addEventListener('DOMContentLoaded', () => {

    // Add click event handler to network parameter tabs.
    document.querySelectorAll('#network-parameter-tab-container button')
            .forEach((tab, index) => {
                const tabId = tab.id;
                const windowId = tabId.replace('-tab', '-parameters');

                tab.addEventListener('click', () => {
                    switchOnParameterWindow(tabId, '#network-parameter-tab-container button', windowId);
                });
            });
    
    // Add click event handler to select-all and reset buttons.
    document.querySelectorAll('.select-all-button')
            .forEach((selectAllButton) => {
                selectAllButton.addEventListener('click', () => {
                    checkTargets(`#${selectAllButton.parentElement.id} input[type=checkbox]`);
                });
            });
            
    document.querySelectorAll('.reset-button')
            .forEach((resetButton) => {
                resetButton.addEventListener('click', () => {
                    uncheckTargets(`#${resetButton.parentElement.id} input[type=checkbox]`);
                });
            });

    // If list_cmpd_classification_superclas is selected, selection of superclass is allowed.
    document.getElementById('id_filter_select_category').addEventListener('change', (e) => {
        const layerFilters = document.querySelectorAll('#layer-filter-container select');
        
        if (e.target.value === 'list_cmpd_classification_superclass') {
            for (layerFilter of layerFilters) {
                if (layerFilter.id === 'chemical-superclass-filer') {
                    layerFilter.style.display = 'inline flow-root';
                    layerFilter.disabled = false;
                } else {
                    layerFilter.style.display = 'none';
                    layerFilter.disabled = true;
                }
            } 
        } else {
            for (layerFilter of layerFilters) {
                layerFilter.disabled = true;
                
                if (layerFilter.id === 'dummy-layer-filter') {
                    layerFilter.style.display = 'inline flow-root';
                } else {
                    layerFilter.style.display = 'none';
                }
            }
        }
    })

    document.getElementById('id_adducts_for_suspect_compound_negative').style.display = 'none';

    // If ion mode is set to positive, only positive adducts will be displayed.
    document.getElementById('id_ion_mode').addEventListener('change', (e) => {
        const posAdductContainer = document.getElementById('id_adducts_for_suspect_compound_positive');
        const negAdductContainer = document.getElementById('id_adducts_for_suspect_compound_negative');
        const posAdducts = posAdductContainer.querySelectorAll('input');
        const negAdducts = negAdductContainer.querySelectorAll('input');

        if (e.target.value === 'pos') {
            posAdductContainer.style.display = '';
            posAdducts.forEach((posAdduct) => {posAdduct.disabled = false});
            negAdductContainer.style.display = 'none';
            negAdducts.forEach((negAdduct) => {negAdduct.disabled = true});
        } else {
            posAdductContainer.style.display = 'none';
            posAdducts.forEach((posAdduct) => {posAdduct.disabled = true});
            negAdductContainer.style.display = '';
            negAdducts.forEach((negAdduct) => {negAdduct.disabled = false});
        }
    });
    
    document.getElementById('id_adducts_for_suspect_compound').setAttribute('expanded', 'false');

    // When #btn-to-expand-adduct-container button is clicked, a list of adducts will be expanded. 
    document.getElementById('btn-to-expand-adduct-container').addEventListener('click', () => {
        const adductContainer = document.getElementById('id_adducts_for_suspect_compound');
        const openIcon = document.querySelector('#btn-to-expand-adduct-container .open-icon');
        const closeIcon = document.querySelector('#btn-to-expand-adduct-container .close-icon');
        
        if (adductContainer.getAttribute('expanded') === 'false') {
            adductContainer.style.height = "max-content";
            adductContainer.setAttribute('expanded', 'true');
            openIcon.style.display = 'none';
            closeIcon.style.display = '';
        } else {
            adductContainer.style.height = "100px";
            adductContainer.setAttribute('expanded', 'false');
            openIcon.style.display = '';
            closeIcon.style.display = 'none';
        }
    });
});

const switchOnParameterWindow = (targetTabId, tabSelectorAll, windowId) => {
    document.querySelectorAll(tabSelectorAll)
            .forEach((tab) => {
                if (tab.id === targetTabId) {
                    tab.classList.add('active');
                } else {
                    tab.classList.remove('active');
                }
            });

    parameterWindows = document.getElementsByClassName('parameter-window');

    for (parameterWindow of parameterWindows) {
        if (parameterWindow.id === windowId) {
            parameterWindow.style.display = '';
        } else {
            parameterWindow.style.display = 'none';
        }
    }
}

const checkTargets = (selectorAll) => {
    document.querySelectorAll(selectorAll)
            .forEach((target) => {
                if (target.type === 'checkbox') {
                    target.checked = true;
                }
            });
}

const uncheckTargets = (selectorAll) => {
    document.querySelectorAll(selectorAll)
            .forEach((target) => {
                if (target.type === 'checkbox') {
                    target.checked = false;
                }
            });
}

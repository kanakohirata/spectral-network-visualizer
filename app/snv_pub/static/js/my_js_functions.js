document.addEventListener('DOMContentLoaded', function(event) {
	console.log('DOMContent is loaded! document.addEventListener() is being called.');

	// Get all form elements.
	const formList = document.getElementsByTagName('form');
	console.log('Number of forms: ' + formList.length);
	// Set validateEmpty() function for submit event on the forms.
	for(let i = 0; i < formList.length; i++) {
    	formList[i].addEventListener('submit', function(e) {
			// Set disabled=true to form controls if their value id empty.
			formTarget = e.target
			for (let i = 0; i < formTarget.elements.length; i++){
				if (! formTarget.elements[i].value){
					formTarget.elements[i].disabled = true;
					console.log(formTarget.elements[i].name + ', disabled=true');
				}
			}
		});
	}
	
	
	// Get elements whose name is 'ion_mode'.
	const ionModeChoices = document.getElementsByName('ion_mode');
	// Set function for click event on the 'ion_mode' elements.
	for(let i = 0; i < ionModeChoices.length; i++) {
    	ionModeChoices[i].addEventListener('click', function() {
			const adductChoicesPos = this.form.elements['prec_type_pos'];
			const adductChoicesNeg = this.form.elements['prec_type_neg'];
			
			// If positive ion mode is selected, enable the positive adduct choices.
			if (this.value === 'pos') {
				adductChoicesPos.style.display = '';
				adductChoicesPos.disabled = false;
				adductChoicesNeg.style.display = 'none';
				adductChoicesNeg.disabled = true;
			} else {
				adductChoicesPos.style.display = 'none';
				adductChoicesPos.disabled = true;
				adductChoicesNeg.style.display = '';
				adductChoicesNeg.disabled = false;
			}
		});
	}
	
	// Get elements whose name is 'ex_mass_prec'.
	const exMassFields = document.getElementsByName('ex_mass_prec');
	// Set function for input event on the 'ex_mass_prec' elements.
	for(let i = 0; i < exMassFields.length; i++) {
    	exMassFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'mass_tol_prec' elements.
			if (this.value) {
				this.form.elements['mass_tol_prec'].required = true;
			} else {
				this.form.elements['mass_tol_prec'].required = false;
			}
		});
	}
	
	// Get elements whose name is 'ex_mass'.
	const exMassCompFields = document.getElementsByName('ex_mass');
	// Set function for input event on the 'ex_mass' elements.
	for(let i = 0; i < exMassCompFields.length; i++) {
		exMassCompFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'mass_tol' elements.
			if (this.value) {
				this.form.elements['mass_tol'].required = true;
			} else {
				this.form.elements['mass_tol'].required = false;
			}
		});
	}
	
	// Get elements whose name is 'ex_mass_frag'.
	const exMassSubFields = document.getElementsByName('ex_mass_frag');
	// Set function for input event on the 'ex_mass_frag' elements.
	for(let i = 0; i < exMassSubFields.length; i++) {
    	exMassSubFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'mass_tol_frag' elements.
			if (this.value) {
				//console.log('mass_tol_frag is required')
				this.form.elements['mass_tol_frag'].required = true;
			} else {
				//console.log('mass_tol_frag is not required')
				this.form.elements['mass_tol_frag'].required = false;
			}
		});
	}
	
	
	// Get elements whose name is 'peaks_mz_vs_ccs'.
	let peaksFields = document.querySelectorAll('[name="peaks_mz_vs_ccs"], [name="peaks_mz_vs_ccs_vs_int"]');
	// Set function for change event on the 'peaks_mz_vs_ccs' elements.
	for(let i = 0; i < peaksFields.length; i++) {
    	peaksFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'mz_tol' and 'ccs_tol' elements.
			if (this.value) {
				this.form.elements['ccs_tol'].required = true;
				return;
			} else {
				this.form.elements['ccs_tol'].required = false;
			}
		});
	}
	
	// Get elements whose name is 'peaks_mz' or 'peaks_mz_vs_ccs'.
	peaksFields = document.querySelectorAll('[name="peaks_mz"], [name="peaks_mz_vs_ccs"], [name="peaks_mz_vs_int"], [name="peaks_mz_vs_ccs_vs_int"]');
	// Set function for change event on the 'peaks_mz' elements.
	for(let i = 0; i < peaksFields.length; i++) {
    	peaksFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'mz_tol' elements.
			if (this.value) {
				this.form.elements['mz_tol'].required = true;
				this.form.elements['cutoff'].required = true;
				// If 'peaks_mz' field is filled, set disabled = true on 'peaks_mz_vs_ccs' field.
				if (peaksFields[i].name === 'peaks_mz') {
					const peaksMzVsCcsFields = document.getElementsByName('peaks_mz_vs_ccs');
					for (let n = 0; n < peaksMzVsCcsFields.length; n++) {
						peaksMzVsCcsFields[n].disabled = true;
					}
				}
				// If 'peaks_mz_vs_ccs' field is filled, set disabled = true on 'peaks_mz' field.
				else if (peaksFields[i].name == 'peaks_mz_vs_ccs') {
					const peaksMzFields = document.getElementsByName('peaks_mz');
					for (let n = 0; n < peaksMzFields.length; n++) {
						peaksMzFields[n].disabled = true;
					}
				}
				// If 'peaks_mz_vs_int' field is filled, set disabled = true on 'peaks_mz_vs_ccs_vs_int' field.
				else if (peaksFields[i].name === 'peaks_mz_vs_int') {
					const peaksMzVsCcsVsIntFields = document.getElementsByName('peaks_mz_vs_ccs_vs_int');
					for (let n = 0; n < peaksMzVsCcsVsIntFields.length; n++) {
						peaksMzVsCcsVsIntFields[n].disabled = true;
					}
				}
				// If 'peaks_mz_vs_ccs_vs_int' field is filled, set disabled = true on 'peaks_mz_vs_int' field.
				else if (peaksFields[i].name === 'peaks_mz_vs_ccs_vs_int') {
					const peaksMzVsCcsVsIntFields = document.getElementsByName('peaks_mz_vs_int');
					for (let n = 0; n < peaksMzVsCcsVsIntFields.length; n++) {
						peaksMzVsCcsVsIntFields[n].disabled = true;
					}
				}
				return;
			} else {
				this.form.elements['mz_tol'].required = false;
				console.log(peaksFields[i].name);
				// If 'peaks_mz' field is blank, set disabled = false on 'peaks_mz_vs_ccs' field.
				if (peaksFields[i].name == 'peaks_mz') {
					const peaksMzVsCcsFields = document.getElementsByName('peaks_mz_vs_ccs');
					for (let n = 0; n < peaksMzVsCcsFields.length; n++) {
						peaksMzVsCcsFields[n].disabled = false;
					}
				}
				// If 'peaks_mz_vs_ccs' field is blank, set disabled = false on 'peaks_mz' field.
				else if (peaksFields[i].name == 'peaks_mz_vs_ccs') {
					const peaksMzFields = document.getElementsByName('peaks_mz');
					for (let n = 0; n < peaksMzFields.length; n++) {
						peaksMzFields[n].disabled = false;
					}
				}
				// If 'peaks_mz_vs_int' field is blank, set disabled = false on 'peaks_mz_vs_ccs_vs_int' field.
				else if (peaksFields[i].name == 'peaks_mz_vs_int') {
					const peaksMzFields = document.getElementsByName('peaks_mz_vs_ccs_vs_int');
					for (let n = 0; n < peaksMzFields.length; n++) {
						peaksMzFields[n].disabled = false;
					}
				}
				// If 'peaks_mz_vs_ccs_vs_int' field is blank, set disabled = false on 'peaks_mz_vs_int' field.
				else if (peaksFields[i].name == 'peaks_mz_vs_ccs_vs_int') {
					const peaksMzFields = document.getElementsByName('peaks_mz_vs_int');
					for (let n = 0; n < peaksMzFields.length; n++) {
						peaksMzFields[n].disabled = false;
					}
				}
			}
		});
	}	
	
	// Get elements whose name is 'prec_mz'.
	const precMzFields = document.getElementsByName('prec_mz');
	// Set function for change event on the 'prec_mz' elements.
	for(let i = 0; i < precMzFields.length; i++) {
    	precMzFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'prec_tol' elements.
			if (this.value) {
				this.form.elements['prec_tol'].required = true;
			} else {
				this.form.elements['prec_tol'].required = false;
			}
		});
	}
	
	// Get elements whose name is 'frag_mz'.
	const fragMzFields = document.getElementsByName('frag_mz');
	// Set function for change event on the 'frag_mz' elements.
	for(let i = 0; i < fragMzFields.length; i++) {
    	fragMzFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'frag_tol' elements.
			if (this.value) {
				this.form.elements['frag_tol'].required = true;
			} else {
				this.form.elements['frag_tol'].required = false;
			}
		});
	}
	
	// Get elements whose name is 'frag_ccs'.
	const fragCcsFields = document.getElementsByName('frag_ccs');
	// Set function for change event on the 'frag_ccs' elements.
	for(let i = 0; i < fragCcsFields.length; i++) {
		fragCcsFields[i].addEventListener('input', function() {
			//If the value is not empty, set required=true to 'frag_ccs_tol' elements.
			if (this.value) {
				this.form.elements['frag_ccs_tol'].required = true;
			} else {
				this.form.elements['frag_ccs_tol'].required = false;
			}
		});
	}
	
	const btnCompSearch = document.getElementById('btn-comp-search');
	if (btnCompSearch) {
		btnCompSearch.click();
	}
	
	const btnFragPeakSearch = document.getElementById('btn-frag-peak-search');
	if (btnFragPeakSearch) {
		btnFragPeakSearch.click();
	}
	
	// Select IMS grah to be displayed in detail page.
	const imsGraphSelect = document.getElementById('ims-graph-select');
	console.log(imsGraphSelect);
	if (imsGraphSelect) {
		if (imsGraphSelect.value == 'mblg') {
			document.getElementById("mobilogram-3d").style.display = "";
			document.getElementById("ims-heatmap").style.display = "none";
		} else if (imsGraphSelect.value == "hmap") {
			document.getElementById("ims-heatmap").style.display = "";
			document.getElementById("mobilogram-3d").style.display = "none";
		}
	}
	
	if (imsGraphSelect) {
		imsGraphSelect.addEventListener('change', (event) => {
			if (imsGraphSelect.value == 'mblg') {
				document.getElementById("mobilogram-3d").style.display = "";
				document.getElementById("ims-heatmap").style.display = "none";
			} else if (imsGraphSelect.value == "hmap") {
				document.getElementById("ims-heatmap").style.display = "";
				document.getElementById("mobilogram-3d").style.display = "none";
			}
		});
	}
	
	// Select a CCS calibration.
	const ccsCalibrationSelect = document.getElementById('ccs-calibration-select');
	console.log(ccsCalibrationSelect);
	if (ccsCalibrationSelect) {
		const ccsCalibrationId = ccsCalibrationSelect.value;
		const ccsElements = document.getElementsByClassName('ccs');
		for (let i = 0; i < ccsElements.length; i++) {
			if (ccsElements[i].classList.contains(`ccs-calibration-${ccsCalibrationId}`)) {
				ccsElements[i].style.display = '';
			} else {
				ccsElements[i].style.display = 'none';
			}
		}
	}
	
	if (ccsCalibrationSelect) {
		ccsCalibrationSelect.addEventListener('change', (event) => {
			const ccsCalibrationId = ccsCalibrationSelect.value;
			const ccsElements = document.getElementsByClassName('ccs');
			for (let i = 0; i < ccsElements.length; i++) {
				if (ccsElements[i].classList.contains(`ccs-calibration-${ccsCalibrationId}`)) {
					ccsElements[i].style.display = '';
				} else {
					ccsElements[i].style.display = 'none';
				}
			}
		});
	}

	// Get elements whose class name contains 'filter-param'.
	const filterCheckBoxes = document.getElementsByClassName('filter-param');
	// console.log('Number of checkbox for filtering: ' + filterCheckBoxes.length);	
	// Set filter parameters in sessionStrage.
	if (filterCheckBoxes.length) {
		let keysOfFilterCheckState = new Set();

		for(let i = 0; i < filterCheckBoxes.length; i++) {
			// console.log(filterCheckBoxes[i]);
			let key = `${filterCheckBoxes[i].id}-checked`;
			keysOfFilterCheckState.add(key)
	
			filterCheckBoxes[i].addEventListener('change', function(e) {
				if (filterCheckBoxes[i].checked) {
					sessionStorage.setItem(key, 1);
					filterCheckBoxes[i].value = 1;
				} else {
					sessionStorage.setItem(key, 0);
					filterCheckBoxes[i].value = 0;
				}
				
				console.log(`${key} = ${sessionStorage.getItem(key)}`);
			});
			
			const checked = parseInt(sessionStorage.getItem(key));
			// console.log(`${key} = ${checked}`);
			if (checked) {
				filterCheckBoxes[i].checked = true;
				filterCheckBoxes[i].value = 1;
			} else {
				filterCheckBoxes[i].checked = false;
				filterCheckBoxes[i].value = 0;
			}
		}
		sessionStorage.setItem('pathnameOfFilterPage', location.pathname);
		sessionStorage.setItem('strKeysOfFilterCheckState', Array.from(keysOfFilterCheckState).join(','));
	}

	// Remove filter parameters of fragments from sessionStorage.
	const pathname = location.pathname;
	if (pathname !== sessionStorage.getItem('pathnameOfFilterPage') & sessionStorage.getItem('strKeysOfFilterCheckState') !== '') {
		const keysOfFilterCheckState = sessionStorage.getItem('strKeysOfFilterCheckState').split(',');
		// console.log('Remove filter parameters of fragments from sessionStorage.')
		for (let i=0; i<keysOfFilterCheckState.length; i++) {
			// console.log(keysOfFilterCheckState[i]);
			sessionStorage.removeItem(keysOfFilterCheckState[i]);
		}
		// console.log(sessionStorage);
	}
});


function switchForms(formId) {
	
	const formSwitches = document.getElementsByClassName('form-switch');
	for (let i = 0; i < formSwitches.length; i++){
		formSwitches[i].classList.remove('active')
	}
	var obj = event.target;
	obj.classList.add('active')
	
	const targetFormClass = 'search-form'
	console.log(targetFormClass);
	formList = document.getElementsByClassName(targetFormClass);
	for (let i = 0; i < formList.length; i++){
		if (formList[i].id === formId){
			formList[i].style.display = '';
			for (let n = 0; n < formList[i].elements.length; n++){
				formList[i].elements[n].disabled = false;
			}
			for (let n = 0; n < formList[i].elements.length; n++){
				if (formList[i].elements[n].defaultChecked == true) {
					formList[i].elements[n].click();
				}
			}
			
		} else {
			formList[i].style.display = 'none';
			for (let n = 0; n < formList[i].elements.length; n++){
				formList[i].elements[n].disabled = true;
			}
		}
	}
}

window.onpageshow = function(event) {
    if (
           event.persisted
        || window.performance && window.performance.navigation.type == 2
    ) {
        window.location.reload();
    }
};

function removeElem(elemId) {
		    const elem = document.getElementById(elemId);
			if (elem) {
				elem.remove();
			}
		}

function wait100Ms() {
	return new Promise((resolve) => {
		setTimeout(() => {
			resolve();
		}, 100);
	});
}

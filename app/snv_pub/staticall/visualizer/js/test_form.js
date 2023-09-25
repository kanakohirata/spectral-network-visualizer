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

const csrfToken = getCookie('csrftoken');

document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('spectral-network-form').addEventListener('submit', async (e) => {
        e.preventDefault();
        let body = new FormData();
        const url = e.target.action;
        const inputParameters = document.querySelectorAll('.network-parameter:not(input[name="adducts_for_suspect_compound"])');

        for (const parameter of inputParameters) {
            if (parameter.value !== '') {
                body.append(parameter.name, parameter.value);
            }
        }

        fetch(url, {
            method: 'POST',
            body: body,
            headers: {
                'X-CSRFToken': csrfToken
            },
        }).then((response) => {
            console.log(response);
            if (!response.ok) {
                throw new Error(response.statusText);
            };
            return response
        }).then((response) => {
            const content_type = response.headers.get('Content-Type');
            console.log(content_type);
            if (content_type.includes('text/html')) {
                response.text().then((text) => {
                    const dom = new DOMParser().parseFromString(text, "text/html");
                    console.log(dom);

                    // document.getElementById('main-content-container').remove();
                    let body = document.getElementById('body');
                    body.replaceWith(dom.activeElement);
                });
                
            } else if (content_type.includes('application/json')) {
                console.log('Get a JsonResponse');
            }

            
        })
    });
});
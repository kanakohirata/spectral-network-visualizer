document.addEventListener('DOMContentLoaded', () => {
    document.querySelectorAll('.confirm-before-submit')
        .forEach((element) => {
            element.addEventListener('click', (event) => {
                let message = element.dataset.confirm;
                if (!message) {
                    message = "Are you sure you want to proceed?";
                }
                if (window.confirm(message)) {
                    ;
                } else {
                    event.preventDefault();
                }
            });
        });
});


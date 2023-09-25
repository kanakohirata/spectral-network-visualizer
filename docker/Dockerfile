# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3.7-slim

RUN apt-get update && apt-get install -y libxrender1 libxext6

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

# Install 

# Install pip requirements
COPY requirements.txt .
RUN python -m pip install -r requirements.txt

WORKDIR /app
# COPY ../app /app

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

# WORKDIR /app/snv-dev

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
# This command isn't correct, can't find the right path to wsgi
# RUN chmod 744 ./startup.sh
# ENTRYPOINT ["./snv-dev/startup.sh"]
# CMD ["bash", "startup.sh"]

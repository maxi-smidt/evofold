FROM python:3.13-slim
LABEL authors="maximiliansmidt"

WORKDIR /app/backend

COPY requirements.txt ./

RUN pip install --no-cache-dir --upgrade -r requirements.txt

COPY . ./

EXPOSE 8080

CMD ["fastapi", "run", "main.py", "--port", "8080", "--proxy-headers"]
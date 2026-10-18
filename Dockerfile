FROM python:3.12-slim

WORKDIR /study
COPY requirements.lock.txt .
RUN pip install --no-cache-dir -r requirements.lock.txt
COPY . .
RUN pip install --no-deps -e .
CMD ["python", "-m", "pytest", "-q"]

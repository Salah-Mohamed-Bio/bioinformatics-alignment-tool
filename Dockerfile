FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY alignment_logic.py .
COPY api.py .

EXPOSE 5000

CMD ["gunicorn", "--bind", "0.0.0.0:5000", "api:app"]
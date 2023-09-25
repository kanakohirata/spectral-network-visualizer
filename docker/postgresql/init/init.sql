CREATE DATABASE snv_pub WITH TEMPLATE = template0 OWNER = postgres ENCODING = 'UTF8' LC_COLLATE = 'C' LC_CTYPE = 'C';

CREATE USER user1 WITH PASSWORD 'user1';

ALTER ROLE user1 SET client_encoding TO 'utf8';

ALTER ROLE user1 SET default_transaction_isolation TO 'read committed';

ALTER ROLE user1 SET timezone TO 'Asia/Tokyo';

ALTER USER user1 CREATEDB;

GRANT ALL PRIVILEGES ON DATABASE snv_pub TO user1;

\c snv_pub;

CREATE SCHEMA django;

GRANT ALL PRIVILEGES ON SCHEMA django TO user1;

\q
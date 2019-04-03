# 1. install rdkit-postgresql
- [ref](https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment)

```bash
conda install -c rdkit rdkit-postgresql
[conda folder]/envs/my-rdkit-env/bin/initdb -D db
pg_ctl -D db/ -l logfile start
```
# 2. config
- [ref](https://www.rdkit.org/docs/Cartridge.html#configuration)
```bash
synchronous_commit = off        # immediate fsync at commit
full_page_writes = off          # recover from partial page writes
shared_buffers = 2048MB         # min 128kB
                                # (change requires restart)
work_mem = 128MB                # min 64kB
```

# 3. create db

# 4. merge db
- [ref](https://stackoverflow.com/questions/9499227/postgresql-merge-2-similar-databases)
```bash
    pg_dump -d db1 -t table1 |psql db2

    then psql and do

    insert into table2 (select * from table1);
```
## dump with gzip [gist](https://gist.github.com/brock/7a7a70300096632cec30)
```
pg_dump -d raw --no-owner -Z 9 > raw.sql.gz
gunzip < raw.sql.gz | psql zinc15
```

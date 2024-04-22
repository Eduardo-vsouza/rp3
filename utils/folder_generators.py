import os


class Content:
    def __init__(self, file, fullfile, group, db, main_dir):
        self.mainDir = main_dir

        self.group = group
        self.groupDir = f'{main_dir}/{group}'

        self.db = db
        self.dbDir = f'{self.groupDir}/{db}'

        self.file = file
        self.fullFile = fullfile


def group_folder_generator(main_dir):
    groups = os.listdir(main_dir)
    for group in groups:
        group_dir = f'{main_dir}/{group}'
        if os.path.isdir(group_dir):

            databases = os.listdir(group_dir)
            for db in databases:
                if db.endswith("target_database.fasta") or db.endswith("_target_decoy.fasta") or db.endswith("target_decoy_database.fasta") or db == 'db':
                    if os.path.isdir(f'{group_dir}/{db}'):
                        db_dir = f'{group_dir}/{db}'
                        files = os.listdir(db_dir)
                        for file in files:
                            fullfile = f'{db_dir}/{file}'
                            content = Content(file=file, fullfile=fullfile, db=db, group=group, main_dir=main_dir)
                            yield content

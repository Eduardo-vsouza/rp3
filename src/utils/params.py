import os


class Params:
    def __init__(self, args):
        self.outdir = args.outdir

        self.modes = ['ribocov', 'database', 'search', 'postms', 'quant']

        self.totalParams = {}
        self.paramsFile = f'{self.outdir}/params.txt'

        self.mode = None

        self.__define_params_file()

    def add_mode_parameters(self, class_instance, args):
        self.mode = class_instance.mode
        params = class_instance.params
        if self.mode not in self.totalParams:
            self.totalParams[self.mode] = params
        else:
            self.totalParams[self.mode].append(params)
        self.totalParams[self.mode].append(args)


    def __define_params_file(self):
        if not os.path.exists(self.paramsFile):
            with open(self.paramsFile, 'w') as outfile:
                lines = []
                for mode in self.modes:
                    lines.append(f'# {mode} mode\n')
                outfile.writelines(lines)

    def update_params_file(self):
        with open(self.paramsFile, 'r') as handler:
            lines = handler.readlines()
            i = self.__locate_mode_section(lines, self.mode)
            for param in self.totalParams[self.mode]:
                i += 1
                lines.insert(i, f'{param}\n')
        with open(self.paramsFile, 'w') as outfile:
            outfile.writelines(lines)

    @staticmethod
    def __locate_mode_section(lines, mode):
        ind = None
        for i, line in enumerate(lines):
            if line.startswith(f"# {mode} mode"):
                ind = i
                break
        return ind




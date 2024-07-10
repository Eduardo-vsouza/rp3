import os

import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk

from ..pipeline_config import PipelineStructure


# # Load the available smORF names and corresponding PNG file paths
# # In practice, you would have a more sophisticated way to load your smORF data
# smorf_data = {
#     'smorf1': 'path/to/smorf1.png',
#     'smorf2': 'path/to/smorf2.png'
# }


class ORFGatherer(PipelineStructure):
    def __init__(self, args):
        super().__init__(args=args)

        self.canonicalData = {}
        self.smorfData = {}

    def collect_data(self):
        files = os.listdir(self.contextFiguresDir)
        for file in files:
            if file.endswith(".png"):
                splat = file.split("_context")
                smorf = splat[0]
                genes = splat[1].split("_")
                for gene in genes:
                    self.canonicalData[gene.replace(".png", "")] = f'{self.contextFiguresDir}/{file}'
                self.smorfData[smorf] = f'{self.contextFiguresDir}/{file}'


class PGViz(PipelineStructure):
    def __init__(self, args, root):
        super().__init__(args=args)
        self.root = root
        self.root.title("RpViz")
        self.smorfData, self.orfData = self.__collect_orfs()

        self.create_widgets()

    def __collect_orfs(self):
        gatherer = ORFGatherer(args=self.args)
        gatherer.collect_data()
        smorfs, orfs = gatherer.smorfData, gatherer.canonicalData
        return smorfs, orfs

    def create_widgets(self):
        # Create a combobox for selecting smORFs
        self.smorf_label = tk.Label(self.root, text="Select smORF:")
        self.smorf_label.pack()

        self.smorf_combobox = ttk.Combobox(self.root, values=list(self.smorfData.keys()))
        self.smorf_combobox.pack()
        self.smorf_combobox.bind("<<ComboboxSelected>>", self.on_combobox_select)

        # Create a combobox for selecting canonical ORFs
        self.orf_label = tk.Label(self.root, text="Select canonical ORF:")
        self.orf_label.pack()

        self.orf_combobox = ttk.Combobox(self.root, values=list(self.orfData.keys()))
        self.orf_combobox.pack()
        self.orf_combobox.bind("<<ComboboxSelected>>", self.on_combobox_orf_select)

        self.image_label = tk.Label(self.root)
        self.image_label.pack(fill=tk.BOTH, expand=1)

    def on_combobox_select(self, event):
        smorf_name = self.smorf_combobox.get()
        self.display_image(smorf_name, self.smorfData)

    def on_combobox_orf_select(self, event):
        orf_name = self.orf_combobox.get()
        self.display_image(orf_name, self.orfData)

    # def create_widgets(self):
    #     # Create a combobox for selecting smORFs
    #     self.smorf_label = tk.Label(self.root, text="Select smORF:")
    #     self.smorf_label.pack()
    #
    #     self.smorf_combobox = ttk.Combobox(self.root, values=list(self.smorfData.keys()))
    #     self.smorf_combobox.pack()
    #     self.smorf_combobox.bind("<<ComboboxSelected>>", self.on_combobox_select)
    #     self.__create_canonical_button()
    #
    #     self.image_label = tk.Label(self.root)
    #     self.image_label.pack(fill=tk.BOTH, expand=1)
    #
    #
    # def __create_canonical_button(self):
    #     self.orfLabel = tk.Label(self.root, text="Select canonical ORF:")
    #     self.orfLabel.pack()
    #
    #     self.orfComboBox = ttk.Combobox(self.root, values=list(self.orfData.keys()))
    #     self.orfComboBox.pack()
    #     self.orfComboBox.bind("<<ComboboxSelected>>", self.on_combobox_orf_select)
    #
    # def on_combobox_select(self, event):
    #     smorf_name = self.smorf_combobox.get()
    #     self.display_image(smorf_name)
    #
    # def on_combobox_orf_select(self, event):
    #     orf_name = self.orfComboBox.get()
    #     self.display_image(orf_name)

    # def display_image(self, smorf_name):
    #     # Retrieve the image file path from smorf_data
    #     image_path = self.smorfData.get(smorf_name)
    #     if image_path:
    #         img = Image.open(image_path)
    #         # img = img.resize((500, 400), Image.ANTIALIAS)  # Resize the image if needed
    #         img_tk = ImageTk.PhotoImage(img)
    #
    #         # Update the image label with the new image
    #         self.image_label.configure(image=img_tk)
    #         self.image_label.image = img_tk
    #     else:
    #         tk.messagebox.showerror("Error", f"No image found for smORF: {smorf_name}")
    def display_image(self, name, data):
        # Retrieve the image file path from the provided data
        image_path = data.get(name)
        if image_path:
            img = Image.open(image_path)
            img_tk = ImageTk.PhotoImage(img)

            # Update the image label with the new image
            self.image_label.configure(image=img_tk)
            self.image_label.image = img_tk
        else:
            messagebox.showerror("Error", f"No image found for: {name}")

# if __name__ == "__main__":
#     root = tk.Tk()
#     interface = PGViz(root)
#     root.mainloop()
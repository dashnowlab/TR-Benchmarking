import gzip
import io


class Reader:
    def __init__(self, file_path: str, buffer_size = io.DEFAULT_BUFFER_SIZE):
        """
        Docstring for __init__
        
        :param file_path: Description
        :type file_path: str
        :param buffer_size: Description
        """
        self.file_obj = None
        self.buffer = buffer_size
        self.path = file_path
        self._raw_line = None
        self.cur_line = None # will be the same as raw_line if no format function is provided to read()
        self.cur_loc = None


    @property
    def raw_line(self):
        """
        Docstring for raw_line
        
        
        """
        return self._raw_line


    def close_file(self):
        """
        Docstring for close_file
        
        
        """
        if not self.file_obj: # if file was opened
            self.file_obj.close()

        self.cur_loc = None


    def open_file(self):
        """
        Docstring for open_file
        
        
        """
        try:
            if self.path.endswith(".gz"):
                self.file_obj = gzip.open(self.path, "rt", encoding="utf-8")
            else:
                self.file_obj = open(self.path, "r", encoding="utf-8", buffering=self.buffer)

            self.cur_loc = 0
            return self
        except (IOError, OSError) as e: 
            raise FileIOError(f"File Opening Error: {e}")


    def read(self, format = None):
        """
        Docstring for read
              
        :param format: Description
        """
        try:
            self._raw_line = self.file_obj.readline()
            self.cur_loc = self.file_obj.tell()

            # if the line is not empty (ie. the end of the file has not been reached) 
            if self._raw_line:
                self.prev_line = self.cur_line
                if format is not None:
                    # then format and set the current line
                    self.cur_line = format(self._raw_line)
                else:
                    self.cur_line = self._raw_line
            else: 
                self.prev_line = self.cur_line
                self.cur_line = self._raw_line

            return self.cur_line
                    
        except gzip.BadGzipFile:
            raise FileReadError(f"Failed to read from {self.path}\nInvalid .gz")
        except UnicodeError:
            raise FileReadError(f"Failed to read from {self.path}\nContains Invalid UTF-8 Characters")
        except Exception as e:
            raise FileReadError(f"Failed to read from {self.path}\n Unknown Error {e}")


    def skipMetaData(self, delimiter='#', end_delimiter = None): 
        """
        Docstring for skipMetaData
        
        
        """

        if self._raw_line is None:
            self.read()

        # loop while the current line contains meta data
        while self._raw_line and self._raw_line.startswith(delimiter):
            if end_delimiter and self._raw_line.startswith(end_delimiter):
                break
            else:
                self.read()  

        # save the file position of the end of the header/metadata
        self.header_end = self.file_obj.tell()


    def _setFilePosition(self, file_pos: int):
        self.file_obj.seek(file_pos)


    def __iter__(self):
        return self
    
    def __next__(self):
        line = self.read(self.formatLine())

        if not line:
            return StopIteration

        return line

    def __enter__(self):
        return self.open_file()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close_file()



class VCFReader(Reader):
    _DEFAULT = object()

    def __init__(self, file_path: str):
        """
        Docstring for __init__
        
        
        :param file_path: Description
        :type file_path: str
        """
        super().__init__(file_path)
        self.prev_line = None
        self.header_end = None

        self.chrom = None
        self.pos = None
        self.end_pos = None
        self.id = None
        self.ref = None
        self.alt = None
        self.qual = None
        self.filter = None
        self.info = None
        self.format = None
        self.sample = None
        self.genotype = None
        

    def buildGt(self, sample_col=9, ref=None, alt = None):
        """
        Docstring for buildGt
        
        
        :param sample_col: Description
        :param ref: Description
        :param alt: Description
        """
        try:
            # use class parameters unless otherwise specified
            ref = self.ref if not ref else ref
            alt = self.alt if not alt else alt
            
            # split line into list for easy data grabbing
            ls = self._raw_line.strip().split("\t")

            # grab genotype indices from current line
            sample_str = ls[sample_col]
            idx_str = sample_str.split(':')[0] # splits sample column and grabs genotype index from start
            choices = [self.ref, *self.alt, None] # list of possible choices to use for constructing genotype
            gt = []

            # loop through the index string and construct genotype
            for idx in idx_str[::2]: # grab characters every 2 indices to skip gt seperator
                # check indices to make sure they are ints (eg. accounting for 1|., etc.) before adding allele to gt
                gt_idx = self._checkIdx(idx) 
                gt.append(choices[gt_idx]) 

            self.genotype = gt

            # add to list of genotypes for comparison between VCFs
            return gt  
        
        except (IndexError, ValueError):
            raise VCFFormatError(f"Failed to construct Genotype using '{idx_str}' from sample: {sample_str}")

        except Exception as e:
            raise VCFFormatError(f"Failed to construct Genotype due to unknown error: {e}\nFrom sample: {sample_str}")


    def formatLine(self, ls: str):
        """
        Docstring for formatLine
        
        
        :param ls: Description
        :type ls: str
        """
        # split line string into list of strings
        line_list = ls.strip().split("\t")

        # if the file is not reading the header/metadata
        if not line_list[0].startswith("#"):  
            try:       
                pos = int(line_list[1])
                ref_len = len(line_list[3])                
                info_col = line_list[7].split(';') # grab the INFO column
                end_str = ""
                for i, str in enumerate(info_col): # seach for end position marker      
                    if "END=" in str:
                        end_str = info_col[i].removeprefix("END=")

                if end_str.isdigit():
                    end_pos = int(end_str)
                else:
                    end_pos = pos + ref_len - 1
                
                # set specific data to their own parameters for better accessibility
                self.chrom = line_list[0]             
                self.pos = pos 
                self.end_pos = end_pos
                self.id = line_list[2]
                self.ref = line_list[3] 
                self.alt = line_list[4].split(",") # returns a list of all alt alleles
                self.qual = line_list[5] 
                self.filter = line_list[6] 
                self.info = line_list[7]
                    
            except ValueError:
                raise VCFFormatError(f"Failed to set position '{line_list[4]}' from line: {line_list}\n")
            
            except IndexError:
                raise VCFFormatError(f"Missing parameter data from line: {line_list}")
            
            except Exception as e:
                raise VCFFormatError(f"Unexpected error setting parameters from line: {line_list}\n{e}")

        return line_list
    

    def read(self, format_method = None):
        if format_method is self._DEFAULT:
            target_format = self.formatLine
        else:
            target_format = format_method
            
        return super().read(target_format)
        

    def _checkIdx(self, idx):
        """
        Docstring for _checkIdx
        
        
        :param idx: Description
        """
        if idx == '.':
            return -1 # will set the allele to None
        else:
            return int(idx)

             


class BEDReader(Reader):
    _DEFAULT = object()

    def __init__(self, file_path: str):
        """
        Docstring for __init__
        
        
        :param file_path: Description
        :type file_path: str
        """
        super().__init__(file_path)
        self.prev_line = None

   
    def formatLine(self, ls: str):
        """
        Docstring for formatLine
        
        
        :param ls: Description
        :type ls: str
        """
        # split line string into list
        line_list = ls.strip().split("\t")

        try:
            # set position information parameters for easier access
            self.chrom = line_list[0]             
            self.pos =  int(line_list[1]) 
            self.end_pos = int(line_list[2])
            self.ref = line_list[3]

        except ValueError:
            raise BEDFormatError(f"ERROR: From file: {self.path}\nFailed to set position {line_list[self.fields['POS']]} from line: {line_list}\n")
        except IndexError:
            raise BEDFormatError(f"ERROR: From file: {self.path}\nMissing parameter data from line: {line_list}")
        except Exception as e:
            raise BEDFormatError(f"ERROR:From file: {self.path}\nUnexpected error setting parameters from line:{line_list}: {e}")


        return line_list
    

    def read(self, format_method = _DEFAULT):
        """
        Docstring for read
              
        :param format_method: Description
        """
        if format_method is self._DEFAULT:
            target_format = self.formatLine
        else:
            target_format = format_method
            
        return super().read(target_format)
    



class FileIOError(Exception):
    def __init__(self, message):
        super().__init__(message)

class FileReadError(Exception):
    def __init__(self, message):
        super().__init__(message)

class BEDFormatError(Exception):
    def __init__(self, message):
        super().__init__(message)

class VCFFormatError(Exception):
    def __init__(self, message):
        super().__init__(message)
class ImgReq:
    image_url = None

    @classmethod
    def set_url(cls, url):
        cls.image_url = url # Store image url

    @classmethod
    def get_url(cls):
        return {"image_url": cls.image_url} if cls.image_url else {"error": "No image URL available!"} # Return the image_url